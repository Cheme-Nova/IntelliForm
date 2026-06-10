"""
modules/mobo_optimizer.py
Multi-Objective Bayesian Optimization (qNEHVI) for IntelliForm.

Replaces the single hand-weighted composite objective used by
modules/bayesian_optimizer.py with a true 4-objective Bayesian
optimization loop:

  1. Cost (USD/kg)        — minimize
  2. Bio-based %          — maximize
  3. Performance score    — maximize
  4. EcoScore (modules.ecometrics composite — biodegradability,
     carbon footprint, ecotoxicity, renewability, REACH/regulatory)
                           — maximize

EcoScore is computed entirely from columns already present in the
ingredient database (no live regulatory API calls), so it can sit
inside the optimization inner loop.

Backend: BoTorch qNEHVI (q-Noisy Expected Hypervolume Improvement) —
the current standard for expensive, noisy, multi-objective Bayesian
optimization. MIT licence. Builds on GPyTorch/PyTorch, already pulled
in by chemprop/gauche.

Reference:
  Daulton et al., "Differentiable Expected Hypervolume Improvement for
  Parallel Multi-Objective Bayesian Optimization", NeurIPS 2020/2021.
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

from modules.bayesian_optimizer import (
    _blend_to_features,
    _composite_objective,
    sample_feasible_candidates_with_relaxation,
)
from modules.ecometrics import compute_ecometrics

try:
    import torch
    from botorch.fit import fit_gpytorch_mll
    from botorch.models import ModelListGP, SingleTaskGP
    from botorch.models.transforms import Normalize, Standardize
    from botorch.acquisition.multi_objective import qNoisyExpectedHypervolumeImprovement
    from botorch.utils.multi_objective.hypervolume import Hypervolume
    from botorch.utils.multi_objective.pareto import is_non_dominated
    from gpytorch.mlls import SumMarginalLogLikelihood
    MOBO_OK = True
except Exception:
    MOBO_OK = False


# ── Data schemas ────────────────────────────────────────────────────────────

@dataclass
class MOBOResult:
    success:              bool
    blend:                Dict[str, float]
    cost_per_kg:          float
    bio_pct:              float
    perf_score:           float
    eco_score:            float
    hypervolume:          float
    n_observations:       int
    acquisition_function: str             = "qNEHVI"
    next_suggestion:      Optional[Dict]  = None
    error_msg:            Optional[str]   = None


@dataclass
class MOBOState:
    """Persistent state — grows with each formulation run."""
    X_observed:    List[List[float]] = field(default_factory=list)  # tabular features
    Y_observed:    List[List[float]] = field(default_factory=list)  # [-cost, bio, perf, eco]
    blend_history: List[Dict]        = field(default_factory=list)
    n_iterations:  int               = 0


FEATURE_COLS_DEFAULT = [
    "Cost_USD_kg", "Bio_based_pct", "Performance_Score",
    "Biodegradability", "Ecotoxicity_Score", "Renewability_Score",
]


# ── Objective vector ─────────────────────────────────────────────────────────

def _objective_vector(blend: Dict[str, float], db: pd.DataFrame,
                       cost: float, bio: float, perf: float) -> Tuple[float, List[float]]:
    """Returns (eco_score, [-cost, bio, perf, eco_score]) — all framed for maximization."""
    eco = compute_ecometrics(blend, db)
    eco_score = eco.eco_score if eco is not None else 50.0
    return eco_score, [-cost, bio, perf, eco_score]


def _reference_point(max_cost: float, min_bio: float, min_perf: float) -> "torch.Tensor":
    """
    Fixed reference point for hypervolume, anchored to the problem's
    constraint floor/ceiling: [-max_cost, min_bio, min_perf, 0].

    A fixed reference point (vs. re-inferring it each iteration) keeps
    hypervolume values comparable across iterations — i.e. non-decreasing
    as the optimizer finds better blends — and gives qNEHVI a stable
    target instead of a moving one.
    """
    eps = 1e-3
    return torch.tensor(
        [-max_cost - eps, min_bio - eps, min_perf - eps, 0.0 - eps],
        dtype=torch.double,
    )


def _hypervolume(Y: np.ndarray, ref_point: "torch.Tensor") -> float:
    """Hypervolume of the non-dominated front of Y vs. a fixed reference point."""
    if not MOBO_OK or len(Y) == 0:
        return 0.0
    train_y = torch.tensor(Y, dtype=torch.double)
    pareto_mask = is_non_dominated(train_y)
    pareto_y = train_y[pareto_mask]
    if pareto_y.shape[0] == 0:
        return 0.0
    hv = Hypervolume(ref_point=ref_point)
    return float(hv.compute(pareto_y))


# ── Main entry point ──────────────────────────────────────────────────────────

def run_mobo_optimization(
    db:            pd.DataFrame,
    max_cost:      float,
    min_bio:       float,
    min_perf:      float,
    state:         Optional[MOBOState] = None,
    n_random_init: int   = 10,
    n_candidates:  int   = 200,
    max_conc:      float = 0.70,
    vertical:      str   = "all",
) -> Tuple[MOBOResult, MOBOState]:
    """
    Run one iteration of multi-objective Bayesian optimization (qNEHVI).

    Phase 1 (< n_random_init iterations): random exploration — spreads
    initial observations across the feasible space for GP fitting.

    Phase 2 (>= n_random_init): qNEHVI-guided exploitation across four
    objectives (cost, bio%, performance, EcoScore).

    Returns (MOBOResult, updated MOBOState).
    """
    if not MOBO_OK:
        return MOBOResult(
            success=False, blend={}, cost_per_kg=0, bio_pct=0, perf_score=0,
            eco_score=0, hypervolume=0, n_observations=0,
            error_msg="BoTorch not available — install botorch.",
        ), state or MOBOState()

    if state is None:
        state = MOBOState()

    idx   = db.set_index("Ingredient") if "Ingredient" in db.columns else db
    names = db["Ingredient"].tolist()

    if len(names) < 3:
        return MOBOResult(
            success=False, blend={}, cost_per_kg=0, bio_pct=0, perf_score=0,
            eco_score=0, hypervolume=0, n_observations=state.n_iterations,
            error_msg="Not enough ingredients in filtered database.",
        ), state

    feature_cols = [c for c in FEATURE_COLS_DEFAULT if c in db.columns]

    candidates = sample_feasible_candidates_with_relaxation(
        db, idx, names, max_cost, min_bio, min_perf, max_conc, n_candidates)

    if not candidates:
        return MOBOResult(
            success=False, blend={}, cost_per_kg=0, bio_pct=0, perf_score=0,
            eco_score=0, hypervolume=0, n_observations=state.n_iterations,
            error_msg="No feasible candidates found. Loosen constraints.",
        ), state

    ref_point = _reference_point(max_cost, min_bio, min_perf)

    # ── Phase 1: Random exploration ───────────────────────────────────────────
    if state.n_iterations < n_random_init or len(state.Y_observed) < 3:
        blend, cost, bio, perf = candidates[np.random.randint(len(candidates))]
        eco_score, y = _objective_vector(blend, db, cost, bio, perf)

        state.X_observed.append(_blend_to_features(blend, db, feature_cols).tolist())
        state.Y_observed.append(y)
        state.blend_history.append(blend)
        state.n_iterations += 1

        return MOBOResult(
            success=True, blend=blend, cost_per_kg=round(cost, 2),
            bio_pct=round(bio, 1), perf_score=round(perf, 1),
            eco_score=round(eco_score, 1),
            hypervolume=round(_hypervolume(np.array(state.Y_observed), ref_point), 4),
            n_observations=state.n_iterations,
            acquisition_function=(f"Random Exploration "
                                   f"(iter {state.n_iterations}/{n_random_init})"),
        ), state

    # ── Phase 2: qNEHVI-guided exploitation ──────────────────────────────────
    train_X = torch.tensor(np.array(state.X_observed), dtype=torch.double)
    train_Y = torch.tensor(np.array(state.Y_observed), dtype=torch.double)

    cand_features = np.array([
        _blend_to_features(b, db, feature_cols) for b, *_ in candidates])
    cand_X = torch.tensor(cand_features, dtype=torch.double)

    try:
        d = train_X.shape[-1]
        models = []
        for i in range(train_Y.shape[-1]):
            models.append(SingleTaskGP(
                train_X, train_Y[:, i:i + 1],
                input_transform=Normalize(d=d),
                outcome_transform=Standardize(m=1),
            ))
        model = ModelListGP(*models)
        mll = SumMarginalLogLikelihood(model.likelihood, model)
        fit_gpytorch_mll(mll)

        acqf = qNoisyExpectedHypervolumeImprovement(
            model=model, ref_point=ref_point, X_baseline=train_X,
            prune_baseline=True,
        )

        with torch.no_grad():
            acq_vals = acqf(cand_X.unsqueeze(1))  # q=1 batches

        best_idx = int(torch.argmax(acq_vals).item())

    except Exception:
        # Fallback: pick by composite objective (mirrors bayesian_optimizer.py)
        obj_vals = [_composite_objective(b, db) for b, *_ in candidates]
        best_idx = int(np.argmax(obj_vals))

    best_blend, best_cost, best_bio, best_perf = candidates[best_idx]
    eco_score, y = _objective_vector(best_blend, db, best_cost, best_bio, best_perf)

    state.X_observed.append(_blend_to_features(best_blend, db, feature_cols).tolist())
    state.Y_observed.append(y)
    state.blend_history.append(best_blend)
    state.n_iterations += 1

    hv = _hypervolume(np.array(state.Y_observed), ref_point)

    # Next experiment suggestion: a different feasible candidate, for diversity
    next_idx = (best_idx + 1) % len(candidates) if len(candidates) > 1 else None
    next_blend = candidates[next_idx][0] if next_idx is not None else None

    return MOBOResult(
        success=True,
        blend=best_blend,
        cost_per_kg=round(best_cost, 2),
        bio_pct=round(best_bio, 1),
        perf_score=round(best_perf, 1),
        eco_score=round(eco_score, 1),
        hypervolume=round(hv, 4),
        n_observations=state.n_iterations,
        acquisition_function="qNEHVI",
        next_suggestion=next_blend,
    ), state
