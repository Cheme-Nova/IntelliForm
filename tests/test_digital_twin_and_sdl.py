"""
tests/test_digital_twin_and_sdl.py
Unit tests for the Manufacturing Digital Twin and Self-Driving Lab (SDL)
integration modules.

Run: pytest tests/ -v
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import pandas as pd
import pytest


@pytest.fixture
def sample_db():
    return pd.DataFrame({
        'Ingredient':       ['Coco-Glucoside', 'Decyl Glucoside', 'Glycerol', 'Citric Acid', 'D-Sorbitol'],
        'SMILES':           [
            'CCCCCCCCCCCCOC1OC(CO)C(O)C(O)C1O',
            'CCCCCCCCCCOC1OC(CO)C(O)C(O)C1O',
            'OCC(O)CO', 'OC(=O)CC(O)(CC(=O)O)C(=O)O',
            'OC[C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO'
        ],
        'Cost_USD_kg':      [3.8, 4.5, 1.5, 1.6, 2.2],
        'Bio_based_pct':    [98,  95,  100, 100, 100],
        'Performance_Score':[88,  82,  65,  68,  72],
        'Stock_kg':         [1850, 920, 5200, 3800, 2100],
        'REACH_Flag':       ['Green','Green','Green','Green','Green'],
        'Function':         ['Primary Surfactant','Primary Surfactant','Humectant','pH Adjuster','Humectant'],
        'Biodegradability': [95, 94, 98, 99, 99],
        'CarbonFootprint_kgCO2eq': [1.2, 1.4, 0.8, 0.6, 0.7],
        'Ecotoxicity_Score':[8, 9, 9, 9, 9],
        'Renewability_Score':[96, 94, 99, 99, 99],
    })


@pytest.fixture
def sample_blend():
    return {
        'Coco-Glucoside': 8.0,
        'Decyl Glucoside': 5.0,
        'Glycerol': 5.0,
        'Citric Acid': 0.3,
        'D-Sorbitol': 3.0,
    }


# ── Digital Twin tests ───────────────────────────────────────────────────────

class TestDigitalTwin:
    def test_simulate_basic(self, sample_db, sample_blend):
        from modules.digital_twin import simulate_manufacturing
        result = simulate_manufacturing(sample_blend, sample_db, batch_kg=500,
                                          vertical="personal_care")
        assert result.batch_kg == 500
        assert result.scale_tier == "pilot"
        assert result.total_cycle_time_min > 0
        assert result.energy_kwh_per_batch > 0
        assert 0 < result.yield_pct <= 100
        assert len(result.process_steps) > 0

    def test_lab_scale_has_no_risks(self, sample_db, sample_blend):
        from modules.digital_twin import simulate_manufacturing
        result = simulate_manufacturing(sample_blend, sample_db, batch_kg=0.5,
                                          vertical="personal_care")
        assert result.scale_tier == "lab"
        assert result.scale_up_risks == []

    def test_production_scale_has_risks(self, sample_db, sample_blend):
        from modules.digital_twin import simulate_manufacturing
        result = simulate_manufacturing(sample_blend, sample_db, batch_kg=20000,
                                          manufacturing_process="emulsification")
        assert result.scale_tier == "production"
        assert len(result.scale_up_risks) > 0
        assert all(r.risk_level == "High" for r in result.scale_up_risks)

    def test_compare_scales(self, sample_db, sample_blend):
        from modules.digital_twin import compare_scales
        cmp = compare_scales(sample_blend, sample_db, vertical="personal_care")
        assert len(cmp.tiers) == 4
        assert cmp.tiers[0].scale_tier == "lab"
        assert cmp.tiers[-1].scale_tier == "production"
        assert "summary" not in dir(cmp) or isinstance(cmp.summary, str)


# ── SDL integration tests ────────────────────────────────────────────────────

class TestSDLProtocol:
    def test_generate_protocol(self, sample_db, sample_blend):
        from modules.sdl_integration import generate_sdl_protocol
        proto = generate_sdl_protocol(sample_blend, sample_db, batch_scale_g=10,
                                        target_specs={"target_ph_min": 5, "target_ph_max": 6.5})
        assert proto.platform == "opentrons_ot2"
        assert len(proto.transfers) == len(sample_blend)
        total_mass = sum(t.target_mass_g for t in proto.transfers)
        assert abs(total_mass - 10.0) < 0.01
        assert "pH (probe, 25°C)" in proto.measurement_requests
        assert proto.opentrons_snippet is not None
        assert "protocol_api" in proto.opentrons_snippet

    def test_generic_platform_has_no_opentrons_snippet(self, sample_db, sample_blend):
        from modules.sdl_integration import generate_sdl_protocol
        proto = generate_sdl_protocol(sample_blend, sample_db, platform="generic_liquid_handler")
        assert proto.opentrons_snippet is None


class TestSDLClosedLoop:
    def test_ingest_pass(self, sample_db, sample_blend):
        from modules.sdl_integration import ingest_sdl_results
        result = ingest_sdl_results(
            sample_blend, sample_db,
            measured={"measured_ph": 6.0},
            target_specs={"target_ph_min": 5, "target_ph_max": 6.5},
        )
        assert result.status == "pass"
        assert result.failure_detected is None
        assert result.next_blend == sample_blend

    def test_ingest_failure_triggers_reformulation(self, sample_db, sample_blend):
        from modules.sdl_integration import ingest_sdl_results
        result = ingest_sdl_results(
            sample_blend, sample_db,
            measured={"measured_ph": 8.5},
            target_specs={"target_ph_min": 5, "target_ph_max": 6.5},
        )
        assert result.status == "reformulate"
        assert result.failure_detected == "ph_too_high"
        assert result.reformulation is not None
        assert abs(sum(result.next_blend.values()) - 100.0) < 0.01

    def test_run_closed_loop_converges(self, sample_db, sample_blend):
        from modules.sdl_integration import run_closed_loop

        def measurement_fn(blend, iteration):
            # pH starts too high and improves with each iteration
            return {"measured_ph": 8.5 - iteration * 1.0}

        history = run_closed_loop(
            sample_blend, sample_db,
            target_specs={"target_ph_min": 5, "target_ph_max": 6.5},
            measurement_fn=measurement_fn,
            max_iterations=5,
        )
        assert len(history) <= 5
        assert history[-1].status == "pass"
