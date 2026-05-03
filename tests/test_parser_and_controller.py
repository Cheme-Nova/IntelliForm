import os
import sys
import asyncio

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


def test_regex_provider_is_direct_and_classifies_common_edge_prompts(monkeypatch):
    import modules.llm_parser as parser_module

    os.environ["LLM_PROVIDER"] = "regex"
    monkeypatch.setattr(parser_module, "_parse_with_groq", lambda text: (_ for _ in ()).throw(AssertionError("groq called")))
    monkeypatch.setattr(parser_module, "_parse_with_anthropic", lambda text: (_ for _ in ()).throw(AssertionError("anthropic called")))
    monkeypatch.setattr(parser_module, "_parse_with_openai", lambda text: (_ for _ in ()).throw(AssertionError("openai called")))
    monkeypatch.setattr(parser_module, "_parse_with_ollama", lambda text: (_ for _ in ()).throw(AssertionError("ollama called")))

    body_wash = parser_module.parse_request("Sulfate-free clear body wash with mild foam under $5/kg.")
    tablet_coating = parser_module.parse_request("Immediate-release tablet film coating system, water-based.")

    assert body_wash.application_type == "personal_care"
    assert tablet_coating.application_type == "pharmaceutical"


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


def test_controller_marks_constraint_violations_and_skips_reports(monkeypatch):
    import api.controller as controller_module
    from modules.llm_parser import ParseResult
    from modules.optimizer import OptResult

    sample_db = pd.DataFrame(
        {
            "Ingredient": ["A", "B"],
            "Vertical": ["fabric_laundry", "fabric_laundry"],
            "Cost_USD_kg": [1.0, 2.0],
            "Bio_based_pct": [70.0, 72.0],
            "Performance_Score": [75.0, 78.0],
            "Function": ["Fabric Builder", "Fabric Primary Surfactant"],
            "REACH_Flag": ["Green", "Green"],
            "Stock_kg": [100, 100],
            "SMILES": ["O", "CCO"],
            "Biodegradability": [90, 88],
            "CarbonFootprint_kgCO2eq": [1.0, 1.2],
            "Ecotoxicity_Score": [9, 8],
            "Renewability_Score": [70, 72],
        }
    )

    monkeypatch.setattr(controller_module, "load_db", lambda: sample_db)
    monkeypatch.setattr(
        controller_module,
        "parse_request",
        lambda text: ParseResult(
            max_cost=3.5,
            min_bio=85.0,
            min_perf=82.0,
            application_type="fabric_laundry",
            reasoning="Laundry brief detected.",
            parser_backend="regex",
            raw_input=text,
        ),
    )
    monkeypatch.setattr(
        controller_module,
        "run_optimization",
        lambda *args, **kwargs: OptResult(
            success=True,
            blend={"A": 60.0, "B": 40.0},
            cost_per_kg=1.4,
            bio_pct=71.0,
            perf_score=76.2,
            status="Optimal",
            vertical="fabric_laundry",
        ),
    )

    def should_not_run(*args, **kwargs):
        raise AssertionError("downstream report should not run for invalid results")

    monkeypatch.setattr(controller_module, "compute_ecometrics", should_not_run)
    monkeypatch.setattr(controller_module, "get_blend_report", should_not_run)
    monkeypatch.setattr(controller_module, "generate_vertical_regulatory_report", should_not_run)
    monkeypatch.setattr(controller_module, "predict_stability", should_not_run)
    monkeypatch.setattr(controller_module, "calculate_carbon_credits", should_not_run)
    monkeypatch.setattr(controller_module, "run_certification_oracle", should_not_run)
    monkeypatch.setattr(controller_module, "run_agent_swarm", lambda result, parsed: [])

    response = controller_module.controller.run(
        input_text="Concentrated liquid detergent, at least 85% bio-based.",
        vertical="fabric_laundry",
        batch_size=1000,
        opt_mode="auto",
        constraints={},
    )

    assert response["result"]["success"] is False
    assert response["result"]["status"] == "ConstraintViolation"
    assert "bio-based 71.0% below 85.0%" in response["result"]["error_msg"]
    assert response["eco"] is None
    assert response["cert"] is None


def test_controller_skips_reports_for_infeasible_result(monkeypatch):
    import api.controller as controller_module
    from modules.llm_parser import ParseResult
    from modules.optimizer import OptResult

    sample_db = pd.DataFrame(
        {
            "Ingredient": ["A"],
            "Vertical": ["industrial"],
            "Cost_USD_kg": [1.0],
            "Bio_based_pct": [80.0],
            "Performance_Score": [85.0],
            "Function": ["Industrial Solvent"],
        }
    )

    monkeypatch.setattr(controller_module, "load_db", lambda: sample_db)
    monkeypatch.setattr(
        controller_module,
        "parse_request",
        lambda text: ParseResult(5.0, 80.0, 85.0, "industrial", "", "regex", text),
    )
    monkeypatch.setattr(
        controller_module,
        "run_optimization",
        lambda *args, **kwargs: OptResult(False, {}, 0, 0, 0, "Infeasible", error_msg="No feasible blend."),
    )
    monkeypatch.setattr(controller_module, "compute_ecometrics", lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError("eco called")))
    monkeypatch.setattr(controller_module, "run_agent_swarm", lambda result, parsed: [])

    response = controller_module.controller.run("Industrial cleaner", "industrial", 1000, "auto", {})

    assert response["result"]["success"] is False
    assert response["eco"] is None
    assert response["carbon"] is None


def test_public_api_helpers_return_expected_shapes():
    from api.main import optimize_pareto, predict_qsar, reformulate
    from api.models import ParetoRequest, QSARRequest, ReformulateRequest

    async def run_checks():
        pareto = await optimize_pareto(
            ParetoRequest(
                vertical="agricultural",
                constraints={"max_cost": 8, "min_bio": 70, "min_perf": 75},
                n_solutions=6,
            )
        )
        assert "frontier" in pareto
        assert "recommended" in pareto
        assert len(pareto["frontier"]) > 0

        qsar = await predict_qsar(
            QSARRequest(
                smiles=["CCO", "O=C(O)c1ccccc1"],
                properties=["biodegradability", "ecotox", "performance"],
            )
        )
        assert "predictions" in qsar
        assert len(qsar["predictions"]) == 2
        assert "biodegradability" in qsar["predictions"][0]

        reform = await reformulate(
            ReformulateRequest(
                blend={"SLS": 12, "Cocamide DEA": 3, "Water": 85},
                failure_type="viscosity",
                vertical="personal_care",
                constraints={},
            )
        )
        assert "root_cause" in reform
        assert "best_suggestion" in reform

    asyncio.run(run_checks())
