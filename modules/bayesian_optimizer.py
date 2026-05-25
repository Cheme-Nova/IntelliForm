"""
modules/bayesian_optimizer.py
Bayesian Optimization for IntelliForm v1.5

GP surrogate + Expected Improvement / UCB acquisition.

Kernel tier (best available at runtime):
  1. GAUCHE TanimotoKernel — fingerprint-space GP; kernel "understands"
     molecular structure. Blend encoded as weighted-average Morgan FP.
     pip install gauche gpytorch
     MIT · Griffiths et al., NeurIPS 2022 (gauche, leojklarner/gauche)

  2. sklearn Matern GP — tabular weighted-average feature vector fallback.
     No extra deps; works on Streamlit Cloud free tier.

References:
  - Griffiths et al., "GAUCHE", NeurIPS 2022
  - Shahriari et al., "Taking the Human Out of the Loop", IEEE 2016
  - ProcessOptimizer: github.com/novonordisk-research/ProcessOptimizer
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple
import warnings
warnings.filterwarnings("ignore")

# ── RDKit — needed for fingerprint encoding ───────────────────────────────────
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    RDKIT_OK = True
except ImportError:
    RDKIT_OK = False

# ── sklearn GP (tier 2 fallback) ─────────────────────────────────────────────
try:
    from sklearn.gaussian_process import GaussianProcessRegressor
    from sklearn.gaussian_process.kernels import Matern, ConstantKernel, WhiteKernel
    from sklearn.preprocessing import StandardScaler
    from scipy.stats import norm
    SKLEARN_GP_OK = True
except ImportError:
    SKLEARN_GP_OK = False

# ── GAUCHE + GPyTorch (tier 1) ────────────────────────────────────────────────
GAUCHE_OK = False
try:
    import torch
    import gpytorch
    from gpytorch.models import ExactGP
    from gpytorch.likelihoods import GaussianLikelihood
    from gpytorch.means import ConstantMean
    from gpytorch.kernels import ScaleKernel
    from gpytorch.mlls import ExactMarginalLogLikelihood
    from gpytorch.distributions import MultivariateNormal
    from gauche.kernels.fingerprint_kernels import TanimotoKernel
    GAUCHE_OK = True
except Exception:
    pass

BAYES_OK = SKLEARN_GP_OK or GAUCHE_OK

MORGAN_RADIUS = 2
MORGAN_NBITS  = 512


# ── Data schemas ──────────────────────────────────────────────────────────────

@dataclass
class BayesianResult:
    success:              bool
    blend:                Dict[str, float]
    cost_per_kg:          float
    bio_pct:              float
    perf_score:           float
    expected_improvement: float
    uncertainty:          float
    n_observations:       int
    acquisition_function: str
    used_gauche:          bool                  = False
    next_suggestion:      Optional[Dict]        = None
    error_msg:            Optional[str]         = None


@dataclass
class BayesianState:
    """Persistent state — grows with each formulation run."""
    X_observed:     List[List[float]] = field(default_factory=list)  # tabular features
    X_fps_observed: List[List[float]] = field(default_factory=list)  # Morgan FPs (GAUCHE)
    y_observed:     List[float]       = field(default_factory=list)
    blend_history:  List[Dict]        = field(default_factory=list)
    n_iterations:   int               = 0
    gp_model:       Optional[object]  = None
    scaler_X:       Optional[object]  = None


# ── Blend encoding ────────────────────────────────────────────────────────────

def _blend_to_fingerprint(blend: Dict[str, float], db: pd.DataFrame) -> Optional[np.ndarray]:
    """
    Weighted-average Morgan fingerprint for a blend.
    Returns float32 array of shape (MORGAN_NBITS,) or None if RDKit unavailable.

    The Tanimoto kernel operates directly on these continuous fingerprint
    vectors — higher overlap → higher kernel similarity → GP learns the
    chemical structure of the objective landscape, not just tabular correlates.
    """
    if not RDKIT_OK:
        return None
    idx   = db.set_index("Ingredient") if "Ingredient" in db.columns else db
    total = sum(blend.values()) or 100.0
    fp    = np.zeros(MORGAN_NBITS, dtype=np.float32)

    for ing, pct in blend.items():
        if ing not in idx.index:
            continue
        smiles = str(idx.loc[ing, "SMILES"]) if "SMILES" in idx.columns else ""
        if not smiles or smiles == "nan":
            continue
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        ing_fp = np.array(
            AllChem.GetMorganFingerprintAsBitVect(mol, MORGAN_RADIUS, nBits=MORGAN_NBITS),
            dtype=np.float32,
        )
        fp += (pct / total) * ing_fp

    # Return None if we got nothing (all ingredients missing SMILES)
    return fp if fp.sum() > 0 else None


def _blend_to_features(blend: Dict[str, float], db: pd.DataFrame,
                        feature_cols: List[str]) -> np.ndarray:
    """Weighted-average tabular feature vector (sklearn fallback)."""
    idx   = db.set_index("Ingredient") if "Ingredient" in db.columns else db
    total = sum(blend.values()) or 100.0
    feat  = np.zeros(len(feature_cols))
    for ing, pct in blend.items():
        if ing in idx.index:
            w = pct / total
            for i, col in enumerate(feature_cols):
                try:
                    feat[i] += w * float(idx.loc[ing, col])
                except Exception:
                    pass
    return feat


# ── Objective function ────────────────────────────────────────────────────────

def _composite_objective(blend: Dict[str, float], db: pd.DataFrame,
                          w_cost: float = 0.3, w_bio: float = 0.4,
                          w_perf: float = 0.3) -> float:
    idx   = db.set_index("Ingredient") if "Ingredient" in db.columns else db
    total = sum(blend.values()) or 100.0
    cost = bio = perf = 0.0
    for ing, pct in blend.items():
        if ing in idx.index:
            w = pct / total / 100
            try:
                cost += w * float(idx.loc[ing, "Cost_USD_kg"])
                bio  += w * float(idx.loc[ing, "Bio_based_pct"])
                perf += w * float(idx.loc[ing, "Performance_Score"])
            except Exception:
                pass
    cost_score = max(0, 1 - cost / 20.0)
    return w_cost * cost_score + w_bio * (bio / 100) + w_perf * (perf / 100)


# ── GAUCHE GP (tier 1) ────────────────────────────────────────────────────────

class _TanimotoGP(ExactGP):
    """Exact GP with GAUCHE Tanimoto fingerprint kernel."""
    def __init__(self, train_x, train_y, likelihood):
        super().__init__(train_x, train_y, likelihood)
        self.mean_module  = ConstantMean()
        self.covar_module = ScaleKernel(TanimotoKernel())

    def forward(self, x):
        return MultivariateNormal(self.mean_module(x), self.covar_module(x))


def _fit_gauche_gp(
    X_fps: np.ndarray,
    y: np.ndarray,
    n_iter: int = 100,
    lr: float = 0.05,
) -> Tuple[object, object]:
    """
    Fit a GAUCHE TanimotoKernel GP on blend fingerprints.
    Returns (model, likelihood) — both in eval mode.
    """
    train_x = torch.tensor(X_fps, dtype=torch.float32)
    train_y = torch.tensor(y,     dtype=torch.float32)

    # Normalise targets to zero mean / unit variance for GP stability
    y_mean = train_y.mean()
    y_std  = train_y.std().clamp(min=1e-6)
    train_y_norm = (train_y - y_mean) / y_std

    likelihood = GaussianLikelihood()
    model      = _TanimotoGP(train_x, train_y_norm, likelihood)

    model.train()
    likelihood.train()
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    mll       = ExactMarginalLogLikelihood(likelihood, model)

    for _ in range(n_iter):
        optimizer.zero_grad()
        output = model(train_x)
        loss   = -mll(output, train_y_norm)
        loss.backward()
        optimizer.step()

    model.eval()
    likelihood.eval()

    # Attach normalisation constants so acquisition can un-normalise
    model._y_mean = y_mean
    model._y_std  = y_std
    return model, likelihood


def _gauche_predict(
    X_cand: np.ndarray,
    model,
    likelihood,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    GP posterior mean and std on candidate fingerprints.
    Returns arrays in original (un-normalised) scale.
    """
    test_x = torch.tensor(X_cand, dtype=torch.float32)
    with torch.no_grad(), gpytorch.settings.fast_pred_var():
        pred = likelihood(model(test_x))
    mu  = (pred.mean  * model._y_std + model._y_mean).numpy()
    std = (pred.stddev * model._y_std).numpy()
    return mu, std


# ── sklearn GP (tier 2 fallback) ─────────────────────────────────────────────

def _fit_sklearn_gp(X: np.ndarray, y: np.ndarray) -> Tuple[object, object]:
    scaler   = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    kernel   = (ConstantKernel(1.0, (1e-3, 1e3)) *
                Matern(length_scale=1.0, length_scale_bounds=(1e-2, 1e2), nu=2.5) +
                WhiteKernel(noise_level=0.01, noise_level_bounds=(1e-5, 0.1)))
    gp = GaussianProcessRegressor(
        kernel=kernel, n_restarts_optimizer=5,
        normalize_y=True, random_state=42,
    )
    gp.fit(X_scaled, y)
    return gp, scaler


# ── Acquisition functions ─────────────────────────────────────────────────────

def _ei(mu: np.ndarray, sigma: np.ndarray, y_best: float,
        xi: float = 0.01) -> np.ndarray:
    sigma = np.maximum(sigma, 1e-8)
    z  = (mu - y_best - xi) / sigma
    ei = (mu - y_best - xi) * norm.cdf(z) + sigma * norm.pdf(z)
    ei[sigma < 1e-8] = 0.0
    return ei


def _ucb(mu: np.ndarray, sigma: np.ndarray, kappa: float = 2.0) -> np.ndarray:
    return mu + kappa * sigma


# ── Main entry point ──────────────────────────────────────────────────────────

def run_bayesian_optimization(
    db:            pd.DataFrame,
    max_cost:      float,
    min_bio:       float,
    min_perf:      float,
    state:         Optional[BayesianState] = None,
    n_random_init: int   = 10,
    n_candidates:  int   = 200,
    max_conc:      float = 0.70,
    acquisition:   str   = "EI",
    vertical:      str   = "all",
) -> Tuple[BayesianResult, BayesianState]:
    """
    Run one iteration of Bayesian optimization.

    Phase 1 (< n_random_init iterations): random exploration
    Phase 2 (>= n_random_init): GP-guided exploitation

    Kernel selection (transparent to caller):
      - GAUCHE TanimotoKernel if gauche + gpytorch installed and SMILES available
      - sklearn Matern GP on tabular features otherwise

    Returns (BayesianResult, updated BayesianState).
    """
    if not BAYES_OK:
        return BayesianResult(
            success=False, blend={}, cost_per_kg=0, bio_pct=0, perf_score=0,
            expected_improvement=0, uncertainty=0, n_observations=0,
            acquisition_function=acquisition,
            error_msg="No GP backend available — install scikit-learn or gauche.",
        ), state or BayesianState()

    if state is None:
        state = BayesianState()

    idx   = db.set_index("Ingredient") if "Ingredient" in db.columns else db
    names = db["Ingredient"].tolist()

    if len(names) < 3:
        return BayesianResult(
            success=False, blend={}, cost_per_kg=0, bio_pct=0, perf_score=0,
            expected_improvement=0, uncertainty=0, n_observations=0,
            acquisition_function=acquisition,
            error_msg="Not enough ingredients in filtered database.",
        ), state

    feature_cols = [c for c in
                    ["Cost_USD_kg", "Bio_based_pct", "Performance_Score",
                     "Biodegradability", "Ecotoxicity_Score", "Renewability_Score"]
                    if c in db.columns]

    # ── Generate constraint-feasible candidates ───────────────────────────────
    def _sample_candidates(n: int, cost_relax: float = 1.0,
                           bio_relax: float = 1.0, perf_relax: float = 1.0):
        out = []
        for _ in range(n):
            n_ings   = np.random.randint(2, min(8, len(names) + 1))
            selected = np.random.choice(names, n_ings, replace=False)
            weights  = np.random.dirichlet(np.ones(n_ings))
            weights  = np.minimum(weights, max_conc)
            weights /= weights.sum()
            blend    = {ing: round(w * 100, 1) for ing, w in zip(selected, weights)}
            cost = sum(float(idx.loc[i, "Cost_USD_kg"]) * p / 100
                       for i, p in blend.items() if i in idx.index)
            bio  = sum(float(idx.loc[i, "Bio_based_pct"]) * p / 100
                       for i, p in blend.items() if i in idx.index)
            perf = sum(float(idx.loc[i, "Performance_Score"]) * p / 100
                       for i, p in blend.items() if i in idx.index)
            if (cost <= max_cost * cost_relax and
                    bio  >= min_bio  * bio_relax and
                    perf >= min_perf * perf_relax):
                out.append((blend, cost, bio, perf))
        return out

    candidates = _sample_candidates(n_candidates)
    if not candidates:
        for relax in [(1.1, 0.9, 0.9), (1.2, 0.8, 0.8), (1.4, 0.7, 0.7)]:
            candidates = _sample_candidates(n_candidates, *relax)
            if candidates:
                break

    if not candidates:
        return BayesianResult(
            success=False, blend={}, cost_per_kg=0, bio_pct=0, perf_score=0,
            expected_improvement=0, uncertainty=0,
            n_observations=state.n_iterations, acquisition_function=acquisition,
            error_msg="No feasible candidates found. Loosen constraints.",
        ), state

    # ── Phase 1: Random exploration ───────────────────────────────────────────
    if state.n_iterations < n_random_init or len(state.y_observed) < 3:
        scored = sorted(
            [(b, c, bio, perf, _composite_objective(b, db))
             for b, c, bio, perf in candidates],
            key=lambda x: -x[4],
        )
        blend, cost, bio, perf, obj_val = scored[0]

        state.X_observed.append(
            _blend_to_features(blend, db, feature_cols).tolist())
        fp = _blend_to_fingerprint(blend, db)
        state.X_fps_observed.append((fp.tolist() if fp is not None
                                     else [0.0] * MORGAN_NBITS))
        state.y_observed.append(obj_val)
        state.blend_history.append(blend)
        state.n_iterations += 1

        return BayesianResult(
            success=True, blend=blend, cost_per_kg=round(cost, 2),
            bio_pct=round(bio, 1), perf_score=round(perf, 1),
            expected_improvement=obj_val, uncertainty=0.5,
            n_observations=state.n_iterations, used_gauche=False,
            acquisition_function=(f"Random Exploration "
                                   f"(iter {state.n_iterations}/{n_random_init})"),
        ), state

    # ── Phase 2: GP-guided exploitation ──────────────────────────────────────
    y = np.array(state.y_observed)
    y_best = float(y.max())

    # Pre-compute candidate features
    cand_tabular = np.array([
        _blend_to_features(b, db, feature_cols) for b, *_ in candidates])
    cand_fps     = np.array([
        (_blend_to_fingerprint(b, db) if GAUCHE_OK and RDKIT_OK
         else np.zeros(MORGAN_NBITS, dtype=np.float32))
        or np.zeros(MORGAN_NBITS, dtype=np.float32)
        for b, *_ in candidates
    ])

    used_gauche  = False
    best_idx     = 0
    best_ei_val  = 0.0
    uncertainty  = 0.3
    mu_all       = np.zeros(len(candidates))
    sigma_all    = np.full(len(candidates), 0.3)

    # ── Try GAUCHE Tanimoto GP ────────────────────────────────────────────────
    if GAUCHE_OK and RDKIT_OK and len(state.X_fps_observed) >= 3:
        X_fps_obs = np.array(state.X_fps_observed, dtype=np.float32)
        # Only proceed if observed FPs are non-trivial
        if X_fps_obs.sum() > 0:
            try:
                model, likelihood = _fit_gauche_gp(X_fps_obs, y)
                mu_all, sigma_all = _gauche_predict(cand_fps, model, likelihood)
                state.gp_model = model

                acq = (_ucb(mu_all, sigma_all)
                       if acquisition == "UCB"
                       else _ei(mu_all, sigma_all, y_best))
                best_idx    = int(np.argmax(acq))
                best_ei_val = float(acq[best_idx])
                uncertainty = float(sigma_all[best_idx])
                used_gauche = True
            except Exception:
                used_gauche = False   # fall through to sklearn GP

    # ── sklearn Matern GP fallback ────────────────────────────────────────────
    if not used_gauche and SKLEARN_GP_OK:
        X_tab_obs = np.array(state.X_observed)
        try:
            gp, scaler = _fit_sklearn_gp(X_tab_obs, y)
            state.gp_model = gp
            state.scaler_X = scaler

            X_scaled = scaler.transform(cand_tabular)
            mu_all, sigma_all = gp.predict(X_scaled, return_std=True)

            acq = (_ucb(mu_all, sigma_all)
                   if acquisition == "UCB"
                   else _ei(mu_all, sigma_all, y_best))
            best_idx    = int(np.argmax(acq))
            best_ei_val = float(acq[best_idx])
            uncertainty = float(sigma_all[best_idx])
        except Exception:
            # Last resort: pick by objective score
            obj_vals = [_composite_objective(b, db) for b, *_ in candidates]
            best_idx    = int(np.argmax(obj_vals))
            best_ei_val = float(obj_vals[best_idx])

    best_blend, best_cost, best_bio, best_perf = candidates[best_idx]

    # Update state
    obj_val = _composite_objective(best_blend, db)
    state.X_observed.append(
        _blend_to_features(best_blend, db, feature_cols).tolist())
    fp = _blend_to_fingerprint(best_blend, db)
    state.X_fps_observed.append(
        fp.tolist() if fp is not None else [0.0] * MORGAN_NBITS)
    state.y_observed.append(obj_val)
    state.blend_history.append(best_blend)
    state.n_iterations += 1

    # Next experiment suggestion: highest GP uncertainty
    max_unc_idx = int(np.argmax(sigma_all))
    next_blend  = (candidates[max_unc_idx][0]
                   if max_unc_idx != best_idx else None)

    kernel_label = ("GAUCHE Tanimoto GP" if used_gauche
                    else f"sklearn Matern GP")
    acq_label    = f"{acquisition} ({kernel_label})"

    return BayesianResult(
        success=True,
        blend=best_blend,
        cost_per_kg=round(best_cost, 2),
        bio_pct=round(best_bio, 1),
        perf_score=round(best_perf, 1),
        expected_improvement=round(best_ei_val, 4),
        uncertainty=round(uncertainty, 4),
        n_observations=state.n_iterations,
        acquisition_function=acq_label,
        used_gauche=used_gauche,
        next_suggestion=next_blend,
    ), state
