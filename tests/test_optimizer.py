"""
tests/test_optimizer.py
Unit tests for IntelliForm v0.8 optimizer, EcoMetrics, and Pareto modules.

Run: pytest tests/ -v
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import pandas as pd
import pytest

# ── Fixtures ──────────────────────────────────────────────────────────────────

@pytest.fixture
def sample_db():
    """Minimal 5-ingredient DB for fast tests."""
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


# ── Optimizer tests ───────────────────────────────────────────────────────────

class TestOptimizer:
    def test_feasible_blend(self, sample_db):
        from modules.optimizer import run_optimization
        result = run_optimization(sample_db, max_cost=5.0, min_bio=90.0, min_perf=75.0)
        assert result.success
        assert result.cost_per_kg <= 5.0
        assert result.bio_pct >= 90.0
        assert result.perf_score >= 75.0

    def test_blend_weights_sum_to_100(self, sample_db):
        from modules.optimizer import run_optimization
        result = run_optimization(sample_db, max_cost=5.0, min_bio=90.0, min_perf=75.0)
        assert result.success
        total = sum(result.blend.values())
        assert abs(total - 100.0) < 1.0, f"Blend sums to {total}, expected ~100"

    def test_infeasible_triggers_relaxation(self, sample_db):
        from modules.optimizer import run_optimization
        # Demand impossible: 100% bio and perf > 95 and cost < $1
        result = run_optimization(sample_db, max_cost=1.0, min_bio=100.0, min_perf=95.0)
        # Should either succeed (relaxed) or fail gracefully
        assert isinstance(result.success, bool)
        if not result.success:
            assert result.error_msg is not None
        if result.relaxed:
            assert result.relaxation_rounds > 0

    def test_infeasible_returns_opt_result(self, sample_db):
        from modules.optimizer import run_optimization
        result = run_optimization(sample_db, max_cost=0.1, min_bio=100.0, min_perf=100.0)
        # Should not raise — always returns OptResult
        assert hasattr(result, 'success')
        assert hasattr(result, 'error_msg')


# ── EcoMetrics tests ──────────────────────────────────────────────────────────

class TestEcoMetrics:
    def test_basic_compute(self, sample_db):
        from modules.ecometrics import compute_ecometrics
        blend = {"Coco-Glucoside": 60.0, "Glycerol": 40.0}
        result = compute_ecometrics(blend, sample_db)
        assert result is not None
        assert 0 <= result.eco_score <= 100
        assert result.grade in ("A+", "A", "B", "C", "D")

    def test_all_axes_in_range(self, sample_db):
        from modules.ecometrics import compute_ecometrics
        blend = {"Coco-Glucoside": 50.0, "Decyl Glucoside": 30.0, "Glycerol": 20.0}
        result = compute_ecometrics(blend, sample_db)
        for axis in [result.biodegradability, result.carbon_footprint,
                     result.ecotoxicity, result.renewability, result.regulatory]:
            assert 0 <= axis <= 100, f"Axis out of range: {axis}"

    def test_empty_blend_returns_none(self, sample_db):
        from modules.ecometrics import compute_ecometrics
        result = compute_ecometrics({}, sample_db)
        assert result is None

    def test_vs_baseline_has_all_axes(self, sample_db):
        from modules.ecometrics import compute_ecometrics
        blend = {"Coco-Glucoside": 100.0}
        result = compute_ecometrics(blend, sample_db)
        expected_keys = {"Biodegradability", "Carbon Footprint", "Ecotoxicity", "Renewability", "Regulatory"}
        assert expected_keys == set(result.vs_baseline.keys())

    def test_radar_data_format(self, sample_db):
        from modules.ecometrics import compute_ecometrics, ecometrics_radar_data
        blend = {"Coco-Glucoside": 70.0, "Glycerol": 30.0}
        eco = compute_ecometrics(blend, sample_db)
        radar = ecometrics_radar_data(eco)
        assert "categories" in radar
        assert "intelliform" in radar
        assert "baseline" in radar
        assert len(radar["categories"]) == 5
        assert len(radar["intelliform"]) == 5


# ── Pareto optimizer tests ────────────────────────────────────────────────────

class TestParetoOptimizer:
    def test_weighted_sum_produces_frontier(self, sample_db):
        from modules.pareto_optimizer import _run_weighted_sum
        result = _run_weighted_sum(sample_db, max_cost=5.0, min_bio=90.0, min_perf=75.0)
        assert result.success
        assert len(result.frontier) >= 1
        assert result.recommended is not None

    def test_all_frontier_solutions_feasible(self, sample_db):
        from modules.pareto_optimizer import _run_weighted_sum
        result = _run_weighted_sum(sample_db, max_cost=5.0, min_bio=90.0, min_perf=75.0)
        for sol in result.frontier:
            assert sol.cost_per_kg <= 5.0 + 0.05  # small tolerance
            assert sol.bio_pct >= 90.0 - 0.5
            assert sol.perf_score >= 75.0 - 0.5

    def test_topsis_returns_valid_solution(self, sample_db):
        from modules.pareto_optimizer import run_pareto_optimization
        result = run_pareto_optimization(sample_db, max_cost=5.0, min_bio=90.0, min_perf=75.0)
        assert result.success
        rec = result.recommended
        assert rec is not None
        assert rec.cost_per_kg > 0
        assert rec.bio_pct > 0
        assert len(rec.blend) > 0

    def test_pareto_dataframe_shape(self, sample_db):
        from modules.pareto_optimizer import run_pareto_optimization, pareto_frontier_dataframe
        result = run_pareto_optimization(sample_db, max_cost=5.0, min_bio=90.0, min_perf=75.0)
        df = pareto_frontier_dataframe(result)
        assert not df.empty
        assert "Cost ($/kg)" in df.columns
        assert "Bio-based (%)" in df.columns
        assert "Perf Score" in df.columns

    def test_infeasible_pareto(self, sample_db):
        from modules.pareto_optimizer import run_pareto_optimization
        result = run_pareto_optimization(sample_db, max_cost=0.5, min_bio=100.0, min_perf=99.0)
        assert isinstance(result.success, bool)
        if not result.success:
            assert result.error_msg is not None


# ── CSV data integrity ────────────────────────────────────────────────────────

class TestIngredientDB:
    def test_csv_loads(self):
        df = pd.read_csv(os.path.join(os.path.dirname(__file__), '..', 'data', 'ingredients_db.csv'))
        assert len(df) >= 30, f"Expected ≥30 ingredients, got {len(df)}"

    def test_required_columns(self):
        df = pd.read_csv(os.path.join(os.path.dirname(__file__), '..', 'data', 'ingredients_db.csv'))
        required = ['Ingredient', 'SMILES', 'Cost_USD_kg', 'Bio_based_pct',
                    'Performance_Score', 'Stock_kg', 'REACH_Flag']
        for col in required:
            assert col in df.columns, f"Missing column: {col}"

    def test_numeric_ranges(self):
        df = pd.read_csv(os.path.join(os.path.dirname(__file__), '..', 'data', 'ingredients_db.csv'))
        assert df['Cost_USD_kg'].between(0, 500).all(), "Cost out of range"
        assert df['Bio_based_pct'].between(0, 100).all(), "Bio% out of range"
        assert df['Performance_Score'].between(0, 100).all(), "Performance out of range"
        assert df['Stock_kg'].ge(0).all(), "Negative stock"

    def test_no_empty_smiles(self):
        df = pd.read_csv(os.path.join(os.path.dirname(__file__), '..', 'data', 'ingredients_db.csv'))
        assert df['SMILES'].notna().all(), "Found NaN SMILES"
        assert (df['SMILES'].str.len() > 0).all(), "Found empty SMILES"

    def test_all_reach_flags_valid(self):
        df = pd.read_csv(os.path.join(os.path.dirname(__file__), '..', 'data', 'ingredients_db.csv'))
        valid = {'Green', 'Amber', 'Red'}
        assert df['REACH_Flag'].isin(valid).all(), f"Invalid REACH flags found"
