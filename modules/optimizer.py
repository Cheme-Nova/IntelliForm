"""
modules/optimizer.py
PuLP LP optimizer with constraint relaxation for IntelliForm.
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


def _solve(db, max_cost, min_bio, min_perf):
    prob = pulp.LpProblem("IntelliForm", pulp.LpMinimize)
    names = db['Ingredient'].tolist()
    vars_ = {n: pulp.LpVariable(f"x_{i}", lowBound=0, upBound=1) for i, n in enumerate(names)}
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
        blend = {n: round(pulp.value(v)*100, 1) for n, v in vars_.items() if (pulp.value(v) or 0) > 0.005}
    return status, blend


def _calc_metrics(blend, db):
    idx = db.set_index('Ingredient')
    cost = sum(idx.loc[k,'Cost_USD_kg'] * pct/100 for k, pct in blend.items())
    bio  = sum(idx.loc[k,'Bio_based_pct'] * pct/100 for k, pct in blend.items())
    perf = sum(idx.loc[k,'Performance_Score'] * pct/100 for k, pct in blend.items())
    return round(cost,2), round(bio,1), round(perf,1)


def run_optimization(db, max_cost, min_bio, min_perf) -> OptResult:
    original = (max_cost, min_bio, min_perf)
    relaxation_schedule = [(-3, -5, 1.00), (-5, -5, 1.10)]
    status, blend = _solve(db, max_cost, min_bio, min_perf)
    relaxed, rounds_used = False, 0
    for round_num, (d_bio, d_perf, cost_mult) in enumerate(relaxation_schedule, start=1):
        if status == 'Optimal':
            break
        relaxed, rounds_used = True, round_num
        max_cost = max_cost * cost_mult
        min_bio  = max(min_bio + d_bio, 75.0)
        min_perf = max(min_perf + d_perf, 65.0)
        track("constraints_relaxed", {"round": round_num})
        status, blend = _solve(db, max_cost, min_bio, min_perf)
    if status != 'Optimal':
        track("error_optimization_failed", {"final_status": status})
        return OptResult(success=False, blend={}, cost_per_kg=0, bio_pct=0, perf_score=0,
                         status=status, relaxed=relaxed, relaxation_rounds=rounds_used,
                         error_msg=f"Could not find feasible blend after {rounds_used} relaxation round(s).")
    cost, bio, perf = _calc_metrics(blend, db)
    track("formulation_generated", {"cost_per_kg": cost, "bio_based_pct": bio,
                                     "performance_score": perf, "relaxed": relaxed})
    return OptResult(success=True, blend=blend, cost_per_kg=cost, bio_pct=bio,
                     perf_score=perf, status=status, relaxed=relaxed, relaxation_rounds=rounds_used)
