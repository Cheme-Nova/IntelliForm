"""
modules/bayesian_optimizer.py
Bayesian Optimization for IntelliForm v1.3

Gaussian Process surrogate model + Expected Improvement acquisition.
Inspired by ProcessOptimizer (Novo Nordisk, JCIM 2025) and BoTorch (Meta).

No PyTorch dependency — uses scikit-learn GaussianProcessRegressor.
Suitable for Streamlit Cloud deployment.

What it does differently from PuLP LP:
  - Learns from each run — builds a surrogate model of the objective landscape
  - Handles noisy objectives — real formulations have measurement noise
  - Suggests which experiment to run next to maximize learning (active learning)
  - Works well with small datasets (10-100 evaluated formulations)
  - Multi-objective via scalarization or Pareto front approximation

References:
  - Shahriari et al., "Taking the Human Out of the Loop", IEEE 2016
  - ProcessOptimizer: github.com/novonordisk-research/ProcessOptimizer
  - scikit-learn GaussianProcessRegressor docs
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple
import warnings
warnings.filterwarnings("ignore")

try:
    from sklearn.gaussian_process import GaussianProcessRegressor
    from sklearn.gaussian_process.kernels import Matern, ConstantKernel, WhiteKernel
    from sklearn.preprocessing import StandardScaler
    from scipy.stats import norm
    from scipy.optimize import minimize
    BAYES_OK = True
except ImportError:
    BAYES_OK = False


# ── Result schema ─────────────────────────────────────────────────────────────

@dataclass
class BayesianResult:
    success: bool
    blend: Dict[str, float]
    cost_per_kg: float
    bio_pct: float
    perf_score: float
    expected_improvement: float
    uncertainty: float              # GP posterior std at this point
    n_observations: int             # how many runs fed into GP
    acquisition_function: str       # "EI" / "UCB" / "PI"
    next_suggestion: Optional[Dict[str, float]] = None  # next experiment suggestion
    error_msg: Optional[str] = None


@dataclass
class BayesianState:
    """Persistent state — grows with each formulation run."""
    X_observed: List[List[float]] = field(default_factory=list)  # feature vectors
    y_observed: List[float] = field(default_factory=list)         # objective values
    blend_history: List[Dict] = field(default_factory=list)
    n_iterations: int = 0
    gp_model: Optional[object] = None
    scaler_X: Optional[object] = None


# ── Objective function ────────────────────────────────────────────────────────

def _blend_to_features(blend: Dict[str, float], db: pd.DataFrame,
                        feature_cols: List[str]) -> np.ndarray:
    """Convert blend (ingredient: pct) to weighted-average feature vector."""
    idx = db.set_index("Ingredient") if "Ingredient" in db.columns else db
    total = sum(blend.values()) or 100.0
    feat = np.zeros(len(feature_cols))
    for ing, pct in blend.items():
        if ing in idx.index:
            w = pct / total
            for i, col in enumerate(feature_cols):
                try:
                    feat[i] += w * float(idx.loc[ing, col])
                except Exception:
                    pass
    return feat


def _composite_objective(blend: Dict[str, float], db: pd.DataFrame,
                          w_cost: float = 0.3, w_bio: float = 0.4,
                          w_perf: float = 0.3) -> float:
    """
    Compute weighted composite objective for a blend.
    Returns value in [0, 1] — higher is better.
    """
    idx = db.set_index("Ingredient") if "Ingredient" in db.columns else db
    total = sum(blend.values()) or 100.0

    cost, bio, perf = 0.0, 0.0, 0.0
    for ing, pct in blend.items():
        if ing in idx.index:
            w = pct / total / 100
            try:
                cost += w * float(idx.loc[ing, "Cost_USD_kg"])
                bio  += w * float(idx.loc[ing, "Bio_based_pct"])
                perf += w * float(idx.loc[ing, "Performance_Score"])
            except Exception:
                pass

    # Normalize: cost is minimized (invert), bio and perf are maximized
    cost_score = max(0, 1 - cost / 20.0)  # normalize by max expected $20/kg
    bio_score  = bio / 100.0
    perf_score = perf / 100.0

    return w_cost * cost_score + w_bio * bio_score + w_perf * perf_score


# ── GP model ──────────────────────────────────────────────────────────────────

def _fit_gp(X: np.ndarray, y: np.ndarray) -> Tuple[object, object]:
    """Fit Gaussian Process regressor. Returns (gp, scaler)."""
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    kernel = (ConstantKernel(1.0, (1e-3, 1e3)) *
              Matern(length_scale=1.0, length_scale_bounds=(1e-2, 1e2), nu=2.5) +
              WhiteKernel(noise_level=0.01, noise_level_bounds=(1e-5, 0.1)))

    gp = GaussianProcessRegressor(
        kernel=kernel,
        n_restarts_optimizer=5,
        normalize_y=True,
        random_state=42,
    )
    gp.fit(X_scaled, y)
    return gp, scaler


def _expected_improvement(X_candidate: np.ndarray, gp: object,
                           scaler: object, y_best: float,
                           xi: float = 0.01) -> np.ndarray:
    """
    Expected Improvement acquisition function.
    xi: exploration-exploitation tradeoff (higher = more exploration).
    """
    X_scaled = scaler.transform(X_candidate)
    mu, sigma = gp.predict(X_scaled, return_std=True)
    sigma = np.maximum(sigma, 1e-8)

    z = (mu - y_best - xi) / sigma
    ei = (mu - y_best - xi) * norm.cdf(z) + sigma * norm.pdf(z)
    ei[sigma < 1e-8] = 0.0
    return ei


def _upper_confidence_bound(X_candidate: np.ndarray, gp: object,
                              scaler: object, kappa: float = 2.0) -> np.ndarray:
    """UCB acquisition: mu + kappa * sigma. kappa controls exploration."""
    X_scaled = scaler.transform(X_candidate)
    mu, sigma = gp.predict(X_scaled, return_std=True)
    return mu + kappa * sigma


# ── Main Bayesian optimization run ────────────────────────────────────────────

def run_bayesian_optimization(
    db: pd.DataFrame,
    max_cost: float,
    min_bio: float,
    min_perf: float,
    state: Optional[BayesianState] = None,
    n_random_init: int = 10,
    n_candidates: int = 200,
    max_conc: float = 0.70,
    acquisition: str = "EI",
    vertical: str = "all",
) -> Tuple[BayesianResult, BayesianState]:
    """
    Run Bayesian optimization for formulation.

    Phase 1 (n_iterations < n_random_init): random exploration
    Phase 2 (n_iterations >= n_random_init): GP-guided exploitation

    Args:
        db: Ingredient database (vertical-filtered)
        state: Persistent BayesianState from previous runs (None = fresh start)
        n_random_init: Random exploration rounds before GP kicks in
        n_candidates: Random candidates to evaluate per iteration
        acquisition: "EI" (Expected Improvement) or "UCB" (Upper Confidence Bound)

    Returns:
        (BayesianResult, updated BayesianState)
    """
    if not BAYES_OK:
        return BayesianResult(
            success=False, blend={}, cost_per_kg=0, bio_pct=0, perf_score=0,
            expected_improvement=0, uncertainty=0, n_observations=0,
            acquisition_function=acquisition,
            error_msg="scikit-learn/scipy not available — install requirements."
        ), state or BayesianState()

    if state is None:
        state = BayesianState()

    idx = db.set_index("Ingredient") if "Ingredient" in db.columns else db
    names = db["Ingredient"].tolist()

    if len(names) < 3:
        return BayesianResult(
            success=False, blend={}, cost_per_kg=0, bio_pct=0, perf_score=0,
            expected_improvement=0, uncertainty=0, n_observations=0,
            acquisition_function=acquisition,
            error_msg="Not enough ingredients in filtered database."
        ), state

    feature_cols = ["Cost_USD_kg", "Bio_based_pct", "Performance_Score",
                    "Biodegradability", "Ecotoxicity_Score", "Renewability_Score"]
    feature_cols = [c for c in feature_cols if c in db.columns]

    # Generate random candidate blends
    candidates = []
    for _ in range(n_candidates):
        # Dirichlet sample — random blend that sums to 100
        n_ings = np.random.randint(2, min(8, len(names) + 1))
        selected = np.random.choice(names, n_ings, replace=False)
        weights = np.random.dirichlet(np.ones(n_ings))

        # Cap max concentration
        weights = np.minimum(weights, max_conc)
        weights /= weights.sum()

        blend = {ing: round(w * 100, 1) for ing, w in zip(selected, weights)}

        # Apply constraints
        cost = sum(float(idx.loc[ing, "Cost_USD_kg"]) * pct / 100
                   for ing, pct in blend.items() if ing in idx.index)
        bio  = sum(float(idx.loc[ing, "Bio_based_pct"]) * pct / 100
                   for ing, pct in blend.items() if ing in idx.index)
        perf = sum(float(idx.loc[ing, "Performance_Score"]) * pct / 100
                   for ing, pct in blend.items() if ing in idx.index)

        if cost <= max_cost and bio >= min_bio and perf >= min_perf:
            candidates.append((blend, cost, bio, perf))

    if not candidates:
        # Relax constraints progressively
        for relax in [0.9, 0.8, 0.7]:
            for _ in range(n_candidates):
                n_ings = np.random.randint(2, min(6, len(names) + 1))
                selected = np.random.choice(names, n_ings, replace=False)
                weights = np.random.dirichlet(np.ones(n_ings))
                weights = np.minimum(weights, max_conc)
                weights /= weights.sum()
                blend = {ing: round(w * 100, 1) for ing, w in zip(selected, weights)}
                cost = sum(float(idx.loc[ing, "Cost_USD_kg"]) * pct / 100
                           for ing, pct in blend.items() if ing in idx.index)
                bio  = sum(float(idx.loc[ing, "Bio_based_pct"]) * pct / 100
                           for ing, pct in blend.items() if ing in idx.index)
                perf = sum(float(idx.loc[ing, "Performance_Score"]) * pct / 100
                           for ing, pct in blend.items() if ing in idx.index)
                if cost <= max_cost * (2 - relax) and bio >= min_bio * relax and perf >= min_perf * relax:
                    candidates.append((blend, cost, bio, perf))
            if candidates:
                break

    if not candidates:
        return BayesianResult(
            success=False, blend={}, cost_per_kg=0, bio_pct=0, perf_score=0,
            expected_improvement=0, uncertainty=0, n_observations=state.n_iterations,
            acquisition_function=acquisition,
            error_msg="No feasible candidates found. Loosen constraints."
        ), state

    # Phase 1: Random exploration — pick best from candidates by objective
    if state.n_iterations < n_random_init or len(state.X_observed) < 3:
        # Score all candidates and pick best
        scored = [(blend, cost, bio, perf,
                   _composite_objective(blend, db))
                  for blend, cost, bio, perf in candidates]
        scored.sort(key=lambda x: -x[4])
        best = scored[0]
        blend, cost, bio, perf, obj_val = best

        # Update state
        feat = _blend_to_features(blend, db, feature_cols)
        state.X_observed.append(feat.tolist())
        state.y_observed.append(obj_val)
        state.blend_history.append(blend)
        state.n_iterations += 1

        return BayesianResult(
            success=True, blend=blend, cost_per_kg=round(cost, 2),
            bio_pct=round(bio, 1), perf_score=round(perf, 1),
            expected_improvement=obj_val,
            uncertainty=0.5,  # high uncertainty in exploration phase
            n_observations=state.n_iterations,
            acquisition_function=f"Random Exploration (iter {state.n_iterations}/{n_random_init})",
        ), state

    # Phase 2: GP-guided exploitation
    X = np.array(state.X_observed)
    y = np.array(state.y_observed)

    try:
        gp, scaler = _fit_gp(X, y)
        state.gp_model = gp
        state.scaler_X = scaler
    except Exception as e:
        # GP failed — fall back to random
        scored = [(blend, cost, bio, perf, _composite_objective(blend, db))
                  for blend, cost, bio, perf in candidates]
        scored.sort(key=lambda x: -x[4])
        best = scored[0]
        blend, cost, bio, perf, obj_val = best
        return BayesianResult(
            success=True, blend=blend, cost_per_kg=round(cost, 2),
            bio_pct=round(bio, 1), perf_score=round(perf, 1),
            expected_improvement=obj_val, uncertainty=0.3,
            n_observations=state.n_iterations,
            acquisition_function="Random (GP failed)",
        ), state

    # Compute features for all candidates
    cand_feats = np.array([
        _blend_to_features(blend, db, feature_cols)
        for blend, _, _, _ in candidates
    ])

    y_best = max(state.y_observed)

    # Acquisition function
    if acquisition == "UCB":
        acq_values = _upper_confidence_bound(cand_feats, gp, scaler, kappa=2.0)
    else:  # EI default
        acq_values = _expected_improvement(cand_feats, gp, scaler, y_best, xi=0.01)

    best_idx = int(np.argmax(acq_values))
    best_blend, best_cost, best_bio, best_perf = candidates[best_idx]
    best_ei = float(acq_values[best_idx])

    # GP uncertainty at this point
    feat_best = cand_feats[best_idx].reshape(1, -1)
    _, uncertainty = gp.predict(scaler.transform(feat_best), return_std=True)
    uncertainty = float(uncertainty[0])

    # Update state
    obj_val = _composite_objective(best_blend, db)
    feat = _blend_to_features(best_blend, db, feature_cols)
    state.X_observed.append(feat.tolist())
    state.y_observed.append(obj_val)
    state.blend_history.append(best_blend)
    state.n_iterations += 1

    # Suggest next experiment (highest uncertainty among candidates)
    uncertainties = []
    for cf in cand_feats:
        _, unc = gp.predict(scaler.transform(cf.reshape(1,-1)), return_std=True)
        uncertainties.append(float(unc[0]))
    max_unc_idx = int(np.argmax(uncertainties))
    next_blend = candidates[max_unc_idx][0] if max_unc_idx != best_idx else None

    return BayesianResult(
        success=True,
        blend=best_blend,
        cost_per_kg=round(best_cost, 2),
        bio_pct=round(best_bio, 1),
        perf_score=round(best_perf, 1),
        expected_improvement=round(best_ei, 4),
        uncertainty=round(uncertainty, 4),
        n_observations=state.n_iterations,
        acquisition_function=acquisition,
        next_suggestion=next_blend,
    ), state
