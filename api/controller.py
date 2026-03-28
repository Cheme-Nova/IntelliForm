import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pandas as pd
from api.memory import memory

from modules.llm_parser import parse_request
from modules.optimizer import run_optimization
from modules.pareto_optimizer import run_pareto_optimization
from modules.bayesian_optimizer import run_bayesian_optimization
from modules.ecometrics import compute_ecometrics
from modules.qsar import predict_properties
from modules.regulatory import get_blend_report
from modules.vertical_regulatory import generate_vertical_regulatory_report
from modules.stability import predict_stability
from modules.carbon_credits import calculate_carbon_credits
from modules.certification_oracle import run_certification_oracle
from modules.agents import run_agent_swarm
from modules.reformulation_intelligence import run_reformulation_intelligence

DB_PATH = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "data", "ingredients_db.csv"
)

def load_db():
    return pd.read_csv(DB_PATH)

class IntelliFormController:
    def run(self, input_text, vertical, batch_size, opt_mode, constraints):
        db = load_db()
        parsed = parse_request(input_text)

        max_cost = constraints.get("max_cost", 50.0)
        min_bio = constraints.get("min_bio", 0.5)
        min_perf = constraints.get("min_perf", 0.7)

        result = run_optimization(
            db, max_cost, min_bio, min_perf,
            max_concentration=1.0,
            vertical=vertical
        )
        eco = compute_ecometrics(result.blend, db)
        reg = get_blend_report(result.blend)
        vreg = generate_vertical_regulatory_report(result.blend, db, vertical)
        stability = predict_stability(result.blend, db)
        carbon = calculate_carbon_credits(result.blend, db, batch_size)
        cert = run_certification_oracle(result.blend, db, vertical, result.bio_pct)
        agents = run_agent_swarm(result, parsed)

        memory.record("formulation_generated", result.blend, vertical)

        return {
            "parsed": _serialize(parsed),
            "result": _serialize(result),
            "eco": _serialize(eco),
            "reg": _serialize(reg),
            "vreg": _serialize(vreg),
            "stability": _serialize(stability),
            "carbon": _serialize(carbon),
            "cert": _serialize(cert),
            "agents": agents
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
