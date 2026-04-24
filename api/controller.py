import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from api.memory import memory
from modules.agents import run_agent_swarm
from modules.bayesian_optimizer import run_bayesian_optimization
from modules.carbon_credits import calculate_carbon_credits
from modules.certification_oracle import run_certification_oracle
from modules.ecometrics import compute_ecometrics
from modules.llm_parser import parse_request
from modules.optimizer import OptResult, run_optimization
from modules.pareto_optimizer import run_pareto_optimization
from modules.regulatory import get_blend_report
from modules.stability import predict_stability
from modules.vertical_regulatory import generate_vertical_regulatory_report
from modules.verticals import filter_db_by_vertical, get_vertical_constraints


DB_PATH = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "data",
    "ingredients_db.csv",
)


VERTICAL_ALIASES = {
    "pharma": "pharmaceutical",
    "pharmaceutical": "pharmaceutical",
    "personal_care": "personal_care",
    "personal care": "personal_care",
    "cosmetics": "personal_care",
    "home_care": "fabric_laundry",
    "home care": "fabric_laundry",
    "laundry": "fabric_laundry",
    "fabric_laundry": "fabric_laundry",
    "textile": "fabric_laundry",
    "paint": "paint_coatings",
    "coatings": "paint_coatings",
    "paint_coatings": "paint_coatings",
    "food_safe": "food",
    "food": "food",
    "industrial": "industrial",
    "agriculture": "agricultural",
    "agricultural": "agricultural",
    "all": "all",
    "unknown": "unknown",
}


def load_db():
    return pd.read_csv(DB_PATH)


def _canonicalize_vertical(value):
    normalized = str(value or "").strip().lower().replace("-", "_")
    return VERTICAL_ALIASES.get(normalized, normalized or "unknown")


def _normalize_threshold(value, default_value):
    if value is None:
        return default_value
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return default_value
    if 0 < numeric <= 1:
        return round(numeric * 100, 1)
    return round(numeric, 1)


def _merge_constraints(parsed, constraints, vertical):
    base_cost, base_bio, base_perf = get_vertical_constraints(
        vertical,
        parsed.max_cost,
        parsed.min_bio,
        parsed.min_perf,
    )
    constraints = constraints or {}
    max_cost = float(constraints.get("max_cost", base_cost))
    min_bio = _normalize_threshold(constraints.get("min_bio"), base_bio)
    min_perf = _normalize_threshold(constraints.get("min_perf"), base_perf)
    return round(max_cost, 2), round(min_bio, 1), round(min_perf, 1)


def _apply_brief_filters(db, input_text):
    text = str(input_text or "").lower()
    filtered = db.copy()

    def drop_matching(pattern):
        nonlocal filtered
        regex = filtered["Ingredient"].astype(str).str.contains(pattern, case=False, regex=True)
        if "Function" in filtered.columns:
            regex = regex | filtered["Function"].astype(str).str.contains(pattern, case=False, regex=True)
        filtered = filtered.loc[~regex].copy()

    if "silicone-free" in text:
        drop_matching(r"silicone|siloxane|dimethicone|cyclomethicone|cyclopentasiloxane")

    if "phosphate-free" in text:
        drop_matching(r"phosphate|polyphosphate")

    if "boron-free" in text:
        drop_matching(r"boron|borate")

    if "low-voc" in text or "high flash point" in text or "high-flash-point" in text:
        drop_matching(r"ethanol|isopropanol|isopropyl|acetone|mek|mibk|ethyl acetate|butyl acetate")

    if any(term in text for term in ["shampoo", "conditioner", "lotion", "cleanser", "serum", "body lotion", "sensitive skin"]):
        drop_matching(r"\(industrial\)|\(agri\)|\(fabric\)|\(food\)|\(pharma")

    if any(term in text for term in ["conditioner", "lotion", "body lotion"]):
        drop_matching(r"ethanol|isopropanol|isopropyl")

    if any(term in text for term in ["shampoo", "conditioner", "lotion", "cleanser"]) and not any(
        term in text for term in ["aha", "bha", "acid", "exfol", "peel", "ph "]
    ):
        drop_matching(r"glycolic acid|pyruvic acid|lactic acid")

    return filtered


def _opt_result_from_solution(solution, status, vertical):
    return OptResult(
        success=True,
        blend=solution.blend,
        cost_per_kg=solution.cost_per_kg,
        bio_pct=solution.bio_pct,
        perf_score=solution.perf_score,
        status=status,
        vertical=vertical,
    )


def _failed_opt_result(error_msg, status, vertical):
    return OptResult(
        success=False,
        blend={},
        cost_per_kg=0,
        bio_pct=0,
        perf_score=0,
        status=status,
        error_msg=error_msg,
        vertical=vertical,
    )


class IntelliFormController:
    def run(self, input_text, vertical, batch_size, opt_mode, constraints):
        db = load_db()
        parsed = parse_request(input_text)

        requested_vertical = _canonicalize_vertical(vertical)
        inferred_vertical = _canonicalize_vertical(parsed.application_type)
        resolved_vertical = requested_vertical if requested_vertical not in ("", "unknown") else inferred_vertical
        if resolved_vertical == "unknown":
            resolved_vertical = inferred_vertical if inferred_vertical != "unknown" else "all"

        filtered_db = filter_db_by_vertical(db, resolved_vertical)
        filtered_db = _apply_brief_filters(filtered_db, input_text)
        max_cost, min_bio, min_perf = _merge_constraints(parsed, constraints, resolved_vertical)

        optimization_mode = str(opt_mode or "auto").lower()
        pareto = None
        bayesian = None

        if optimization_mode == "pareto":
            pareto = run_pareto_optimization(filtered_db, max_cost, min_bio, min_perf)
            if pareto.success and pareto.recommended:
                result = _opt_result_from_solution(pareto.recommended, pareto.backend, resolved_vertical)
            else:
                result = _failed_opt_result(
                    pareto.error_msg or "Pareto optimization failed.",
                    pareto.backend,
                    resolved_vertical,
                )
        elif optimization_mode == "bayesian":
            bayesian, _ = run_bayesian_optimization(
                filtered_db,
                max_cost,
                min_bio,
                min_perf,
                vertical=resolved_vertical,
            )
            if bayesian.success:
                result = OptResult(
                    success=True,
                    blend=bayesian.blend,
                    cost_per_kg=bayesian.cost_per_kg,
                    bio_pct=bayesian.bio_pct,
                    perf_score=bayesian.perf_score,
                    status="bayesian",
                    vertical=resolved_vertical,
                )
            else:
                result = _failed_opt_result(
                    bayesian.error_msg or "Bayesian optimization failed.",
                    "bayesian",
                    resolved_vertical,
                )
        else:
            result = run_optimization(
                filtered_db,
                max_cost,
                min_bio,
                min_perf,
                max_concentration=1.0,
                vertical=resolved_vertical,
            )

        eco = compute_ecometrics(result.blend, filtered_db)
        reg = get_blend_report(result.blend)
        vreg = generate_vertical_regulatory_report(result.blend, filtered_db, resolved_vertical)
        stability = predict_stability(result.blend, filtered_db)
        carbon = calculate_carbon_credits(result.blend, filtered_db, batch_size)
        cert = run_certification_oracle(result.blend, filtered_db, resolved_vertical, result.bio_pct)
        agents = run_agent_swarm(result, parsed)

        if result.success:
            memory.record("formulation_generated", result.blend, resolved_vertical)

        return {
            "parsed": _serialize(parsed),
            "result": _serialize(result),
            "eco": _serialize(eco),
            "reg": _serialize(reg),
            "vreg": _serialize(vreg),
            "stability": _serialize(stability),
            "carbon": _serialize(carbon),
            "cert": _serialize(cert),
            "agents": agents,
            "pareto": _serialize(pareto),
            "bayesian": _serialize(bayesian),
            "meta": {
                "requested_vertical": requested_vertical,
                "inferred_vertical": inferred_vertical,
                "resolved_vertical": resolved_vertical,
                "optimization_mode_requested": optimization_mode,
                "constraints_used": {
                    "max_cost": max_cost,
                    "min_bio": min_bio,
                    "min_perf": min_perf,
                },
                "ingredient_pool_size": len(filtered_db),
            },
        }


def _serialize(obj):
    """Convert dataclasses/objects to dicts for JSON serialization."""
    if obj is None:
        return None
    if hasattr(obj, "__dict__"):
        return {k: _serialize(v) for k, v in obj.__dict__.items()}
    if isinstance(obj, dict):
        return {k: _serialize(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_serialize(i) for i in obj]
    return obj


controller = IntelliFormController()
