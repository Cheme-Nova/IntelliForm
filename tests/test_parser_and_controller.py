import os
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))


def test_regex_parser_normalizes_vertical_and_constraints():
    from modules.llm_parser import parse_request

    os.environ["LLM_PROVIDER"] = "regex"
    result = parse_request(
        "Silicone-free rinse-off hair conditioner for textured hair, naturally derived, under $6/kg."
    )

    assert result.application_type == "personal_care"
    assert result.max_cost == 6.0
    assert result.min_bio >= 95.0
    assert result.min_perf >= 85.0


def test_controller_uses_parsed_constraints_and_filters_vertical(monkeypatch):
    import api.controller as controller_module
    from modules.llm_parser import ParseResult
    from modules.optimizer import OptResult

    sample_db = pd.DataFrame(
        {
            "Ingredient": ["A", "B", "C"],
            "Vertical": ["fabric_laundry", "fabric_laundry", "industrial"],
            "Cost_USD_kg": [1.0, 2.0, 3.0],
            "Bio_based_pct": [95.0, 90.0, 20.0],
            "Performance_Score": [85.0, 84.0, 70.0],
            "Function": ["Fabric Builder", "Fabric Primary Surfactant", "Industrial Solvent"],
            "REACH_Flag": ["Green", "Green", "Amber"],
            "Stock_kg": [100, 100, 100],
            "SMILES": ["O", "CCO", "CCC"],
            "Biodegradability": [90, 88, 40],
            "CarbonFootprint_kgCO2eq": [1.0, 1.2, 4.0],
            "Ecotoxicity_Score": [9, 8, 3],
            "Renewability_Score": [95, 92, 10],
        }
    )

    captured = {}

    monkeypatch.setattr(controller_module, "load_db", lambda: sample_db)
    monkeypatch.setattr(
        controller_module,
        "parse_request",
        lambda text: ParseResult(
            max_cost=3.5,
            min_bio=92.0,
            min_perf=84.0,
            application_type="fabric_laundry",
            reasoning="Laundry brief detected.",
            parser_backend="regex",
            raw_input=text,
        ),
    )

    def fake_run_optimization(db, max_cost, min_bio, min_perf, max_concentration=1.0, vertical="all"):
        captured["ingredients"] = db["Ingredient"].tolist()
        captured["max_cost"] = max_cost
        captured["min_bio"] = min_bio
        captured["min_perf"] = min_perf
        captured["vertical"] = vertical
        return OptResult(
            success=True,
            blend={"A": 60.0, "B": 40.0},
            cost_per_kg=1.4,
            bio_pct=93.0,
            perf_score=84.5,
            status="Optimal",
            vertical=vertical,
        )

    monkeypatch.setattr(controller_module, "run_optimization", fake_run_optimization)
    monkeypatch.setattr(controller_module, "compute_ecometrics", lambda blend, db: None)
    monkeypatch.setattr(controller_module, "get_blend_report", lambda blend: None)
    monkeypatch.setattr(controller_module, "generate_vertical_regulatory_report", lambda blend, db, vertical: None)
    monkeypatch.setattr(controller_module, "predict_stability", lambda blend, db: None)
    monkeypatch.setattr(controller_module, "calculate_carbon_credits", lambda blend, db, batch_size: None)
    monkeypatch.setattr(controller_module, "run_certification_oracle", lambda blend, db, vertical, bio_pct: None)
    monkeypatch.setattr(controller_module, "run_agent_swarm", lambda result, parsed: [])
    monkeypatch.setattr(controller_module.memory, "record", lambda *args, **kwargs: None)

    response = controller_module.controller.run(
        input_text="Cold-water laundry detergent, phosphate-free, under $3.50/kg.",
        vertical="home_care",
        batch_size=1000,
        opt_mode="auto",
        constraints={},
    )

    assert captured["vertical"] == "fabric_laundry"
    assert captured["ingredients"] == ["A", "B"]
    assert captured["max_cost"] == 3.5
    assert captured["min_bio"] >= 92.0
    assert captured["min_perf"] >= 84.0
    assert response["meta"]["resolved_vertical"] == "fabric_laundry"
