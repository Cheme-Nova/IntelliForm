"""
modules/optimizer.py
PuLP linear programming optimizer for IntelliForm v0.7.

Extracted from app.py so it can be unit-tested independently
and reused by any future API endpoint or batch runner.

Constraint relaxation strategy (auto-triggered on Infeasible):
  Round 1: loosen bio by 3 pts, perf by 5 pts
  Round 2: loosen bio by 5 more pts, perf by 5 more pts, raise cost ceiling 10%
  After that: give up and return failed OptResult
"""
from dataclasses import dataclass, field
from typing import Dict, Optional
import pandas as pd
import pulp

from modules.analytics import track


# ── Result schema ─────────────────────────────────────────────────────────────

@dataclass
class OptResult:
    success: bool
    blend: Dict[str, float]          # {ingredient_name: weight_pct}
    cost_per_kg: float
    bio_pct: float
    perf_score: float
    status: str                      # PuLP status string
    relaxed: bool = False            # True if we had to relax constraints
    relaxation_rounds: int = 0
    error_msg: Optional[str] = None


# ── Core solve ────────────────────────────────────────────────────────────────

def _solve(db: pd.DataFrame, max_cost: float, min_bio: float, min_perf: float) -> tuple:
    """Run one PuLP solve. Returns (status_str, blend_dict)."""
    prob = pulp.LpProblem("IntelliForm_v07", pulp.LpMinimize)
    names = db['Ingredient'].tolist()
    vars_ = {n: pulp.LpVariable(f"x_{i}", lowBound=0, upBound=1)
             for i, n in enumerate(names)}

    idx = db.set_index('Ingredient')

    # Objective: minimise cost
    prob += pulp.lpSum(idx.loc[n, 'Cost_USD_kg'] * vars_[n] for n in names)

    # Constraints
    prob += pulp.lpSum(vars_.values()) == 1                                      # weights sum to 1
    prob += pulp.lpSum(idx.loc[n, 'Cost_USD_kg']    * vars_[n] for n in names) <= max_cost
    prob += pulp.lpSum(idx.loc[n, 'Bio_based_pct']  * vars_[n] for n in names) >= min_bio
    prob += pulp.lpSum(idx.loc[n, 'Performance_Score'] * vars_[n] for n in names) >= min_perf

    prob.solve(pulp.PULP_CBC_CMD(msg=0))
    status = pulp.LpStatus[prob.status]

    blend = {}
    if status == 'Optimal':
        blend = {n: round(pulp.value(v) * 100, 1)
                 for n, v in vars_.items() if (pulp.value(v) or 0) > 0.005}
    return status, blend


def _calc_metrics(blend: dict, db: pd.DataFrame) -> tuple:
    idx = db.set_index('Ingredient')
    cost = sum(idx.loc[k, 'Cost_USD_kg']       * pct / 100 for k, pct in blend.items())
    bio  = sum(idx.loc[k, 'Bio_based_pct']     * pct / 100 for k, pct in blend.items())
    perf = sum(idx.loc[k, 'Performance_Score'] * pct / 100 for k, pct in blend.items())
    return round(cost, 2), round(bio, 1), round(perf, 1)


# ── Public interface ──────────────────────────────────────────────────────────

def run_optimization(
    db: pd.DataFrame,
    max_cost: float,
    min_bio: float,
    min_perf: float
) -> OptResult:
    """
    Run PuLP optimization with automatic constraint relaxation on Infeasible.
    Always returns an OptResult — never raises.
    """
    original = (max_cost, min_bio, min_perf)
    relaxation_schedule = [
        # (delta_bio, delta_perf, cost_multiplier)
        (-3,  -5,  1.00),
        (-5,  -5,  1.10),
    ]

    status, blend = _solve(db, max_cost, min_bio, min_perf)

    relaxed = False
    rounds_used = 0

    for round_num, (d_bio, d_perf, cost_mult) in enumerate(relaxation_schedule, start=1):
        if status == 'Optimal':
            break

        relaxed = True
        rounds_used = round_num
        new_cost = max_cost * cost_mult
        new_bio  = max(min_bio  + d_bio,  75.0)
        new_perf = max(min_perf + d_perf, 65.0)

        track("constraints_relaxed", {
            "round":            round_num,
            "original_max_cost": original[0],
            "original_min_bio":  original[1],
            "original_min_perf": original[2],
            "relaxed_max_cost":  new_cost,
            "relaxed_min_bio":   new_bio,
            "relaxed_min_perf":  new_perf,
            "lp_status":         status
        })

        status, blend = _solve(db, new_cost, new_bio, new_perf)
        max_cost, min_bio, min_perf = new_cost, new_bio, new_perf

    if status != 'Optimal':
        track("error_optimization_failed", {
            "final_status":  status,
            "original_max_cost": original[0],
            "original_min_bio":  original[1],
            "original_min_perf": original[2],
            "relaxation_rounds": rounds_used
        })
        return OptResult(
            success=False, blend={}, cost_per_kg=0, bio_pct=0, perf_score=0,
            status=status, relaxed=relaxed, relaxation_rounds=rounds_used,
            error_msg=f"Could not find feasible blend after {rounds_used} relaxation round(s). "
                      f"Try loosening your constraints manually."
        )

    cost, bio, perf = _calc_metrics(blend, db)

    track("formulation_generated", {
        "cost_per_kg":        cost,
        "bio_based_pct":      bio,
        "performance_score":  perf,
        "num_ingredients":    len(blend),
        "top_ingredient":     max(blend, key=blend.get) if blend else "",
        "relaxed":            relaxed,
        "relaxation_rounds":  rounds_used,
        "final_max_cost":     max_cost,
        "final_min_bio":      min_bio,
        "final_min_perf":     min_perf
    })

    return OptResult(
        success=True, blend=blend, cost_per_kg=cost,
        bio_pct=bio, perf_score=perf, status=status,
        relaxed=relaxed, relaxation_rounds=rounds_used
    )
