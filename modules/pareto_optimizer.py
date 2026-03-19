"""
modules/pareto_optimizer.py
Multi-objective Pareto optimization for IntelliForm v0.8.

Replaces the single-objective PuLP min-cost approach with NSGA-III,
producing a Pareto frontier of non-dominated blends across three axes:

  Objective 1: Minimize cost (USD/kg)
  Objective 2: Maximize bio-based % (negated for minimization)
  Objective 3: Maximize performance score (negated for minimization)

Falls back gracefully to the existing PuLP optimizer if pymoo is not
installed (so the app never breaks in constrained environments).

Usage:
    from modules.pareto_optimizer import run_pareto_optimization, ParetoResult
    pareto = run_pareto_optimization(db, max_cost=5.0, min_bio=90, min_perf=80)
    # pareto.frontier is a list of OptResult-compatible dicts
    # pareto.recommended is the single best balanced blend
"""
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple
import numpy as np
import pandas as pd

from modules.optimizer import run_optimization, OptResult


# ── Result schema ─────────────────────────────────────────────────────────────

@dataclass
class ParetoSolution:
    """One non-dominated solution on the Pareto frontier."""
    blend: Dict[str, float]       # {ingredient: weight%}
    cost_per_kg: float
    bio_pct: float
    perf_score: float
    eco_score: Optional[float]    # filled by EcoMetrics if available
    solution_id: int              # index on frontier, 0 = most balanced


@dataclass
class ParetoResult:
    """Full result from multi-objective optimization."""
    success: bool
    frontier: List[ParetoSolution]     # all Pareto-optimal blends
    recommended: Optional[ParetoSolution]  # best trade-off (TOPSIS selection)
    n_solutions: int
    backend: str                       # "nsga3" | "weighted_sum" | "pulp_fallback"
    error_msg: Optional[str] = None


# ── NSGA-III backend (pymoo) ──────────────────────────────────────────────────

def _run_nsga3(
    db: pd.DataFrame,
    max_cost: float,
    min_bio: float,
    min_perf: float,
    n_gen: int = 200,
    pop_size: int = 100,
) -> Optional[ParetoResult]:
    """
    Run NSGA-III multi-objective optimization via pymoo.
    Returns None if pymoo is not installed.
    """
    try:
        from pymoo.algorithms.moo.nsga3 import NSGA3
        from pymoo.core.problem import Problem
        from pymoo.optimize import minimize
        from pymoo.util.ref_dirs import get_reference_directions
        from pymoo.termination import get_termination
    except ImportError:
        return None

    names = db["Ingredient"].tolist()
    n = len(names)
    idx = db.set_index("Ingredient")

    costs  = np.array([idx.loc[nm, "Cost_USD_kg"]       for nm in names])
    bios   = np.array([idx.loc[nm, "Bio_based_pct"]      for nm in names])
    perfs  = np.array([idx.loc[nm, "Performance_Score"]  for nm in names])

    class FormulationProblem(Problem):
        def __init__(self):
            super().__init__(
                n_var=n,
                n_obj=3,
                n_ieq_constr=3,
                xl=np.zeros(n),
                xu=np.ones(n),
            )

        def _evaluate(self, X, out, *args, **kwargs):
            # Normalize rows to sum = 1
            row_sums = X.sum(axis=1, keepdims=True)
            X_norm = X / (row_sums + 1e-12)

            # Objectives: minimize cost, maximize bio (negate), maximize perf (negate)
            f1 = X_norm @ costs          # minimize cost
            f2 = -(X_norm @ bios)        # maximize bio
            f3 = -(X_norm @ perfs)       # maximize perf

            # Constraints (≤ 0 means feasible) — small tolerance for numerical stability
            TOLS = 0.5   # allow 0.5% slack on bio/perf constraints
            g1 = (X_norm @ costs) - max_cost
            g2 = (min_bio  - TOLS) - (X_norm @ bios)
            g3 = (min_perf - TOLS) - (X_norm @ perfs)

            out["F"] = np.column_stack([f1, f2, f3])
            out["G"] = np.column_stack([g1, g2, g3])

    ref_dirs = get_reference_directions("das-dennis", 3, n_partitions=6)
    algorithm = NSGA3(pop_size=max(len(ref_dirs), pop_size), ref_dirs=ref_dirs)
    termination = get_termination("n_gen", n_gen)

    res = minimize(
        FormulationProblem(),
        algorithm,
        termination,
        seed=42,
        verbose=False,
    )

    if res.X is None or len(res.X) == 0:
        return ParetoResult(
            success=False, frontier=[], recommended=None,
            n_solutions=0, backend="nsga3",
            error_msg="NSGA-III found no feasible solutions — try relaxing constraints."
        )

    # Decode solutions
    solutions = []
    for i, x in enumerate(res.X):
        x_norm = x / (x.sum() + 1e-12)
        blend = {
            names[j]: round(float(x_norm[j]) * 100, 1)
            for j in range(n) if x_norm[j] > 0.01
        }
        c = float(x_norm @ costs)
        b = float(x_norm @ bios)
        p = float(x_norm @ perfs)
        solutions.append(ParetoSolution(
            blend=blend,
            cost_per_kg=round(c, 2),
            bio_pct=round(b, 1),
            perf_score=round(p, 1),
            eco_score=None,
            solution_id=i,
        ))

    recommended = _topsis_select(solutions)

    return ParetoResult(
        success=True,
        frontier=solutions,
        recommended=recommended,
        n_solutions=len(solutions),
        backend="nsga3",
    )


# ── Weighted-sum enumeration fallback ────────────────────────────────────────

def _run_weighted_sum(
    db: pd.DataFrame,
    max_cost: float,
    min_bio: float,
    min_perf: float,
    n_points: int = 12,
) -> ParetoResult:
    """
    Approximate Pareto frontier by sweeping weight combinations over PuLP.
    Works without pymoo. Produces n_points non-dominated blends.
    """
    import pulp

    names = db["Ingredient"].tolist()
    idx = db.set_index("Ingredient")
    costs_ = {nm: idx.loc[nm, "Cost_USD_kg"]       for nm in names}
    bios_  = {nm: idx.loc[nm, "Bio_based_pct"]      for nm in names}
    perfs_ = {nm: idx.loc[nm, "Performance_Score"]  for nm in names}

    # Max values for normalization
    max_c  = max(costs_.values())
    max_b  = max(bios_.values())
    max_p  = max(perfs_.values())

    all_solutions: List[ParetoSolution] = []

    # Generate weight combinations (w_cost, w_bio, w_perf) summing to 1
    weights = []
    steps = int(np.sqrt(n_points)) + 1
    for i in range(steps + 1):
        for j in range(steps + 1 - i):
            k_ = steps - i - j
            weights.append((i / steps, j / steps, k_ / steps))

    seen_costs = set()
    for w_c, w_b, w_p in weights:
        prob = pulp.LpProblem("pareto_ws", pulp.LpMinimize)
        v = {nm: pulp.LpVariable(f"x_{nm}", 0, 1) for nm in names}

        # Scalarized objective
        prob += pulp.lpSum(
            (w_c * costs_[nm] / max_c
             - w_b * bios_[nm]  / max_b
             - w_p * perfs_[nm] / max_p) * v[nm]
            for nm in names
        )
        prob += pulp.lpSum(v.values()) == 1
        prob += pulp.lpSum(costs_[nm] * v[nm] for nm in names) <= max_cost
        prob += pulp.lpSum(bios_[nm]  * v[nm] for nm in names) >= min_bio
        prob += pulp.lpSum(perfs_[nm] * v[nm] for nm in names) >= min_perf
        prob.solve(pulp.PULP_CBC_CMD(msg=0))

        if pulp.LpStatus[prob.status] != "Optimal":
            continue

        blend = {nm: round(pulp.value(v[nm]) * 100, 1) for nm in names if (pulp.value(v[nm]) or 0) > 0.01}
        c = round(sum(costs_[nm] * pct / 100 for nm, pct in blend.items()), 2)
        b = round(sum(bios_[nm]  * pct / 100 for nm, pct in blend.items()), 1)
        p = round(sum(perfs_[nm] * pct / 100 for nm, pct in blend.items()), 1)

        # Deduplicate by cost (rounded to 2 dp)
        key = round(c, 2)
        if key in seen_costs:
            continue
        seen_costs.add(key)

        all_solutions.append(ParetoSolution(
            blend=blend, cost_per_kg=c, bio_pct=b, perf_score=p,
            eco_score=None, solution_id=len(all_solutions)
        ))

    if not all_solutions:
        return ParetoResult(
            success=False, frontier=[], recommended=None,
            n_solutions=0, backend="weighted_sum",
            error_msg="No feasible blends found. Try relaxing constraints."
        )

    # Keep only Pareto-non-dominated solutions
    frontier = _pareto_filter(all_solutions)
    for i, s in enumerate(frontier):
        s.solution_id = i

    recommended = _topsis_select(frontier)
    return ParetoResult(
        success=True, frontier=frontier, recommended=recommended,
        n_solutions=len(frontier), backend="weighted_sum"
    )


# ── Pareto filtering ──────────────────────────────────────────────────────────

def _pareto_filter(solutions: List[ParetoSolution]) -> List[ParetoSolution]:
    """Remove dominated solutions. Minimizing cost, maximizing bio and perf."""
    non_dominated = []
    for s in solutions:
        dominated = False
        for other in solutions:
            if (other.cost_per_kg <= s.cost_per_kg and
                    other.bio_pct >= s.bio_pct and
                    other.perf_score >= s.perf_score and
                    (other.cost_per_kg < s.cost_per_kg or
                     other.bio_pct > s.bio_pct or
                     other.perf_score > s.perf_score)):
                dominated = True
                break
        if not dominated:
            non_dominated.append(s)
    return sorted(non_dominated, key=lambda x: x.cost_per_kg)


# ── TOPSIS selection ──────────────────────────────────────────────────────────

def _topsis_select(solutions: List[ParetoSolution]) -> Optional[ParetoSolution]:
    """
    TOPSIS multi-criteria decision: select the Pareto solution with the
    best balance across cost, bio%, and performance.
    """
    if not solutions:
        return None
    if len(solutions) == 1:
        return solutions[0]

    matrix = np.array([
        [s.cost_per_kg, s.bio_pct, s.perf_score]
        for s in solutions
    ], dtype=float)

    # Normalize
    norms = np.linalg.norm(matrix, axis=0)
    norms[norms == 0] = 1
    normed = matrix / norms

    # Weights: cost 40%, bio 30%, perf 30%
    w = np.array([0.40, 0.30, 0.30])
    weighted = normed * w

    # Ideal best: min cost, max bio, max perf
    ideal_best  = np.array([weighted[:, 0].min(), weighted[:, 1].max(), weighted[:, 2].max()])
    ideal_worst = np.array([weighted[:, 0].max(), weighted[:, 1].min(), weighted[:, 2].min()])

    d_best  = np.linalg.norm(weighted - ideal_best,  axis=1)
    d_worst = np.linalg.norm(weighted - ideal_worst, axis=1)

    scores = d_worst / (d_best + d_worst + 1e-12)
    best_idx = int(np.argmax(scores))
    return solutions[best_idx]


# ── Public interface ──────────────────────────────────────────────────────────

def run_pareto_optimization(
    db: pd.DataFrame,
    max_cost: float,
    min_bio: float,
    min_perf: float,
    n_gen: int = 100,
    pop_size: int = 50,
) -> ParetoResult:
    """
    Run multi-objective Pareto optimization.

    Priority:
      1. NSGA-III via pymoo (if installed)
      2. Weighted-sum enumeration via PuLP (always available)

    Returns ParetoResult with full frontier + TOPSIS-recommended blend.
    """
    # Try NSGA-III first
    result = _run_nsga3(db, max_cost, min_bio, min_perf, n_gen=n_gen, pop_size=pop_size)
    if result is not None:
        return result

    # Fall back to weighted-sum
    return _run_weighted_sum(db, max_cost, min_bio, min_perf)


def pareto_frontier_dataframe(result: ParetoResult) -> pd.DataFrame:
    """
    Convert Pareto frontier to a clean DataFrame for display / export.
    """
    if not result.frontier:
        return pd.DataFrame()
    rows = []
    for s in result.frontier:
        top_ing = max(s.blend, key=s.blend.get) if s.blend else "—"
        rows.append({
            "ID":               s.solution_id,
            "Cost ($/kg)":      s.cost_per_kg,
            "Bio-based (%)":    s.bio_pct,
            "Perf Score":       s.perf_score,
            "EcoScore™":        s.eco_score or "—",
            "Top Ingredient":   top_ing,
            "# Ingredients":    len(s.blend),
        })
    return pd.DataFrame(rows)
