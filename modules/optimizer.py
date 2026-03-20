"""
modules/optimizer.py
PuLP LP optimizer with constraint relaxation for IntelliForm v1.0.

New in v1.0: max_concentration parameter — limits any single ingredient
to a maximum fraction of the blend, forcing more diverse results.
"""
from dataclasses import dataclass
from typing import Dict, Optional
import pandas as pd
import pulp
from modules.analytics import track


@dataclass
class OptResult:
    success: bool
    blend: Dict[str, float]
    cost_per_kg: float
    bio_pct: float
    perf_score: float
    status: str
    relaxed: bool = False
    relaxation_rounds: int = 0
    error_msg: Optional[str] = None


def _solve(db, max_cost, min_bio, min_perf, max_concentration=1.0):
    prob = pulp.LpProblem("IntelliForm", pulp.LpMinimize)
    names = db['Ingredient'].tolist()
    # max_concentration is a fraction (0-1), e.g. 0.70 = max 70% per ingredient
    vars_ = {n: pulp.LpVariable(f"x_{i}", lowBound=0, upBound=max_concentration)
             for i, n in enumerate(names)}
    idx = db.set_index('Ingredient')
    prob += pulp.lpSum(idx.loc[n,'Cost_USD_kg'] * vars_[n] for n in names)
    prob += pulp.lpSum(vars_.values()) == 1
    prob += pulp.lpSum(idx.loc[n,'Cost_USD_kg'] * vars_[n] for n in names) <= max_cost
    prob += pulp.lpSum(idx.loc[n,'Bio_based_pct'] * vars_[n] for n in names) >= min_bio
    prob += pulp.lpSum(idx.loc[n,'Performance_Score'] * vars_[n] for n in names) >= min_perf
    prob.solve(pulp.PULP_CBC_CMD(msg=0))
    status = pulp.LpStatus[prob.status]
    blend = {}
    if status == 'Optimal':
        blend = {n: round(pulp.value(v)*100, 1)
                 for n, v in vars_.items() if (pulp.value(v) or 0) > 0.005}
    return status, blend


def _calc_metrics(blend, db):
    idx = db.set_index('Ingredient')
    cost = sum(idx.loc[k,'Cost_USD_kg'] * pct/100 for k, pct in blend.items())
    bio  = sum(idx.loc[k,'Bio_based_pct'] * pct/100 for k, pct in blend.items())
    perf = sum(idx.loc[k,'Performance_Score'] * pct/100 for k, pct in blend.items())
    return round(cost,2), round(bio,1), round(perf,1)


def run_optimization(db, max_cost, min_bio, min_perf,
                     max_concentration: float = 1.0) -> OptResult:
    """
    Run PuLP optimization with optional concentration cap.

    Args:
        max_concentration: max fraction for any single ingredient (0.3-1.0).
                          0.70 = no ingredient can exceed 70% of blend.
                          1.0 = no constraint (original behavior).
    """
    original = (max_cost, min_bio, min_perf)
    relaxation_schedule = [(-3, -5, 1.00), (-5, -5, 1.10)]

    # If concentration cap is very tight, relax it slightly on failure
    conc = max(max_concentration, 0.30)  # never below 30% — would be infeasible

    status, blend = _solve(db, max_cost, min_bio, min_perf, max_concentration=conc)
    relaxed, rounds_used = False, 0

    for round_num, (d_bio, d_perf, cost_mult) in enumerate(relaxation_schedule, start=1):
        if status == 'Optimal':
            break
        relaxed, rounds_used = True, round_num
        max_cost = max_cost * cost_mult
        min_bio  = max(min_bio + d_bio, 75.0)
        min_perf = max(min_perf + d_perf, 65.0)
        # Also relax concentration cap slightly on failure
        conc = min(conc + 0.10, 1.0)
        track("constraints_relaxed", {"round": round_num})
        status, blend = _solve(db, max_cost, min_bio, min_perf, max_concentration=conc)

    if status != 'Optimal':
        track("error_optimization_failed", {"final_status": status})
        return OptResult(
            success=False, blend={}, cost_per_kg=0, bio_pct=0, perf_score=0,
            status=status, relaxed=relaxed, relaxation_rounds=rounds_used,
            error_msg=f"Could not find feasible blend after {rounds_used} relaxation round(s). "
                      f"Try loosening constraints or increasing max ingredient %."
        )

    cost, bio, perf = _calc_metrics(blend, db)
    track("formulation_generated", {
        "cost_per_kg": cost, "bio_based_pct": bio,
        "performance_score": perf, "relaxed": relaxed,
        "max_concentration": max_concentration,
    })
    return OptResult(
        success=True, blend=blend, cost_per_kg=cost, bio_pct=bio,
        perf_score=perf, status=status, relaxed=relaxed,
        relaxation_rounds=rounds_used
    )
