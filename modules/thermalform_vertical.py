"""
modules/thermalform_vertical.py
ThermalForm(TM) vertical for IntelliForm: PFAS-free single-phase immersion cooling fluids.

Runs IntelliForm's NSGA-III Pareto engine (same pymoo NSGA-III, Das-Dennis reference
directions and sparse seeding as modules/pareto_optimizer.py) with ThermalForm physics as
the objectives, then picks the recommended blend with IntelliForm's TOPSIS selector.

Mapping onto IntelliForm's three axes (so _pareto_filter and _topsis_select apply unchanged):
  cost_per_kg -> facility pumping cost, USD per MW-year (operating cost, minimize)
  bio_pct     -> readily biodegradable share of the blend, % (maximize)
  perf_score  -> chip-side cooling score = 100 x R_benchmark / R_blend (maximize)

Constraints: OCP immersion-fluid spec at worst-case ester permittivity
(Dk <= 2.3, flash >= flash_min, pour < -30 C, viscosity at 50 C < 15 mPa s, no mineral oil),
at most 3 base stocks, each at least 5%.

Physics comes from the ThermalForm package (properties, heat sink, loop, OCP checks).
"""
from __future__ import annotations

import math
import os
import sys
from dataclasses import dataclass
from typing import Dict, List, Optional

import numpy as np

from modules.pareto_optimizer import ParetoSolution, ParetoResult, _pareto_filter, _topsis_select

TF_PATH = os.environ.get("THERMALFORM_PATH", "")
if TF_PATH and TF_PATH not in sys.path:
    sys.path.insert(0, TF_PATH)

from thermalform.properties import load_base_stocks, evaluate          # noqa: E402
from thermalform.physics import HeatSink, forced                        # noqa: E402
from thermalform.loop import loop                                        # noqa: E402
from thermalform.ocp import check, DK_LIMIT                              # noqa: E402

Q_CHIP, T_IN, FLOW = 1000.0, 40.0, 6 / 60000
SINK = HeatSink()
MU50_MAX = 0.015
POUR_UNKNOWN = -31.0     # PAO / isoparaffin pour points unpublished; treated as just passing, flagged


@dataclass
class ThermalFormContext:
    bs: object
    R_bench: float
    pump_usd_bench: float
    eligible: List[int]


def _context(bench=("PAO4", 0.88, "PAO2", 0.12)) -> ThermalFormContext:
    bs = load_base_stocks()
    w = np.zeros(bs.n)
    w[bs.ids.index(bench[0])] = bench[1]
    w[bs.ids.index(bench[2])] = bench[3]
    R = forced(w, bs, FLOW, SINK, Q_CHIP, T_IN)["R_fluid"]
    P = loop(w, bs)["pump_usd_yr"]
    eligible = [i for i, c in enumerate(bs.chemistry) if c != "mineral oil"]   # OCP: avoid mineral oils
    return ThermalFormContext(bs, R, P, eligible)


def score_blend(w: np.ndarray, ctx: ThermalFormContext) -> Dict:
    bs = ctx.bs
    r = forced(w, bs, FLOW, SINK, Q_CHIP, T_IN)
    lp = loop(w, bs)
    oc = check(w, bs)
    p = evaluate(w, bs)
    pour = oc["pour_worst_c"] if not math.isnan(oc["pour_worst_c"]) else POUR_UNKNOWN
    from thermalform.properties import fluid_at
    mu50 = fluid_at(w, bs, 50.0)["mu"]
    return {
        "R": r["R_fluid"], "pump_usd_yr": lp["pump_usd_yr"], "pump_kw": lp["pump_kw"],
        "perf_score": 100.0 * ctx.R_bench / r["R_fluid"],
        "bio_pct": 100.0 * p["biodegradable_frac"],
        "dk_worst": oc["dk_ghz_worst"], "flash_c": oc["flash_c"], "pour_c": pour, "mu50": mu50,
        "dT_vs_bench": (r["R_fluid"] - ctx.R_bench) * Q_CHIP,
        "pump_usd_vs_bench": lp["pump_usd_yr"] - ctx.pump_usd_bench,
    }


def _repair(x: np.ndarray, eligible: List[int], max_stocks=3, min_frac=0.05) -> np.ndarray:
    w = np.zeros_like(x)
    xe = np.clip(x[eligible], 0, None)
    if xe.sum() <= 0:
        xe = np.random.random(len(eligible))
    top = np.argsort(xe)[::-1][:max_stocks]
    vals = xe[top] / xe[top].sum()
    keep = vals >= min_frac
    top, vals = top[keep], vals[keep] / vals[keep].sum()
    for t, v in zip(top, vals):
        w[eligible[t]] = v
    return w


def run_thermalform_pareto(flash_min: float = 150.0, n_gen: int = 120, pop_size: int = 120,
                           seed: int = 42) -> tuple[ParetoResult, List[Dict], ThermalFormContext]:
    from pymoo.algorithms.moo.nsga3 import NSGA3
    from pymoo.core.problem import ElementwiseProblem
    from pymoo.core.repair import Repair
    from pymoo.core.sampling import Sampling
    from pymoo.optimize import minimize
    from pymoo.util.ref_dirs import get_reference_directions

    ctx = _context()
    n = ctx.bs.n

    class TFRepair(Repair):
        def _do(self, problem, X, **kwargs):
            return np.array([_repair(x, ctx.eligible) for x in X])

    class SparseSampling(Sampling):      # IntelliForm sparse seeding, 1-3 active stocks
        def _do(self, problem, n_samples, **kwargs):
            X = np.zeros((n_samples, problem.n_var))
            for i in range(n_samples):
                k = np.random.randint(1, 4)
                active = np.random.choice(ctx.eligible, size=k, replace=False)
                X[i, active] = np.random.random(k)
            return X

    class Problem(ElementwiseProblem):
        def __init__(self):
            super().__init__(n_var=n, n_obj=3, n_ieq_constr=4, xl=np.zeros(n), xu=np.ones(n))

        def _evaluate(self, x, out, *args, **kwargs):
            w = x / (x.sum() + 1e-12)
            s = score_blend(w, ctx)
            out["F"] = [s["pump_usd_yr"], -s["bio_pct"], -s["perf_score"]]
            out["G"] = [s["dk_worst"] - DK_LIMIT, flash_min - s["flash_c"], s["pour_c"] + 30.0, s["mu50"] - MU50_MAX]

    ref_dirs = get_reference_directions("das-dennis", 3, n_partitions=12)
    algo = NSGA3(pop_size=max(len(ref_dirs), pop_size), ref_dirs=ref_dirs,
                 sampling=SparseSampling(), repair=TFRepair(), eliminate_duplicates=True)
    res = minimize(Problem(), algo, ("n_gen", n_gen), seed=seed, verbose=False)
    if res.X is None or len(np.atleast_2d(res.X)) == 0:
        return ParetoResult(False, [], None, 0, "nsga3-thermalform", "no feasible blend"), [], ctx

    sols, details = [], []
    for i, x in enumerate(np.atleast_2d(res.X)):
        w = x / x.sum()
        s = score_blend(w, ctx)
        blend = {ctx.bs.ids[j]: round(100 * w[j], 1) for j in np.argsort(w)[::-1] if w[j] > 1e-6}
        sols.append(ParetoSolution(blend=blend, cost_per_kg=round(s["pump_usd_yr"], 0), bio_pct=round(s["bio_pct"], 1),
                                   perf_score=round(s["perf_score"], 2), eco_score=None, solution_id=i))
        details.append({"blend": blend, **{k: (round(v, 4) if isinstance(v, float) else v) for k, v in s.items()}})
    front = _pareto_filter(sols)
    rec = _topsis_select(front)
    return ParetoResult(True, front, rec, len(front), "nsga3-thermalform"), details, ctx
