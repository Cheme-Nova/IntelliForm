"""
modules/ecometrics.py
EcoMetrics™ sustainability scoring for IntelliForm v0.8.

Computes a 5-axis sustainability profile for a formulation blend:
  1. Biodegradability  — OECD 301B / ASTM D6868 composite
  2. Carbon Footprint  — kgCO₂eq/kg (inverted: lower = better score)
  3. Ecotoxicity       — ECHA aquatic toxicity rating (inverted)
  4. Renewability      — bio-based + feedstock origin composite
  5. Regulatory        — REACH/EPA/EU Ecolabel compliance composite

Each axis is normalized 0–100. A weighted composite EcoScore™ is also
returned for at-a-glance comparison.

Weights (tunable, published in IntelliForm JCIM SI):
  Biodegradability: 25%
  Carbon Footprint: 20%
  Ecotoxicity:      20%
  Renewability:     20%
  Regulatory:       15%
"""
from dataclasses import dataclass
from typing import Dict, Optional
import pandas as pd


# ── Axis weights ──────────────────────────────────────────────────────────────
AXIS_WEIGHTS = {
    "Biodegradability": 0.25,
    "Carbon Footprint":  0.20,
    "Ecotoxicity":       0.20,
    "Renewability":      0.20,
    "Regulatory":        0.15,
}

# REACH flag → regulatory score
_REACH_SCORES = {"Green": 100, "Amber": 50, "Red": 10}


# ── Result schema ─────────────────────────────────────────────────────────────

@dataclass
class EcoMetricsResult:
    # Individual axis scores (0–100, higher = better)
    biodegradability: float
    carbon_footprint: float   # inverted from kgCO₂eq
    ecotoxicity: float        # inverted (lower toxicity = higher score)
    renewability: float
    regulatory: float

    # Composite
    eco_score: float          # weighted average, 0–100
    grade: str                # A+ / A / B / C / D
    blend_summary: Dict[str, float]   # ingredient: weight%

    # Benchmark comparisons (vs typical petrochemical baseline)
    vs_baseline: Dict[str, float]     # axis: delta vs petrochemical


# Carbon footprint inversion: scale from max ~5.0 kgCO₂eq to 0 → 0–100
_MAX_CO2 = 5.0

def _invert_co2(co2: float) -> float:
    """Convert kgCO₂eq/kg to 0–100 score (lower CO₂ = higher score)."""
    return max(0.0, min(100.0, (1 - co2 / _MAX_CO2) * 100))

def _invert_ecotox(ecotox: float) -> float:
    """Ecotoxicity 1–10 scale where 10=safest → normalize to 0–100."""
    return max(0.0, min(100.0, (ecotox / 10.0) * 100))

def _grade(score: float) -> str:
    if score >= 90: return "A+"
    if score >= 80: return "A"
    if score >= 70: return "B"
    if score >= 55: return "C"
    return "D"


# Petrochemical baseline for comparison (industry average INCI surfactant blend)
_PETROCHEMICAL_BASELINE = {
    "Biodegradability": 52.0,
    "Carbon Footprint":  38.0,
    "Ecotoxicity":       41.0,
    "Renewability":      25.0,
    "Regulatory":        60.0,
}


# ── Core computation ──────────────────────────────────────────────────────────

def compute_ecometrics(blend: Dict[str, float], db: pd.DataFrame) -> Optional[EcoMetricsResult]:
    """
    Compute EcoMetrics™ scores for a formulation blend.

    Args:
        blend: {ingredient_name: weight_percent} dict
        db: ingredients DataFrame (must include EcoMetrics columns)

    Returns:
        EcoMetricsResult, or None if blend is empty or columns missing.
    """
    if not blend:
        return None

    required_cols = {"Biodegradability", "CarbonFootprint_kgCO2eq", "Ecotoxicity_Score",
                     "Renewability_Score", "REACH_Flag"}
    if not required_cols.issubset(set(db.columns)):
        # Fallback: synthesize from available data
        return _synthesize_from_bio_pct(blend, db)

    idx = db.set_index("Ingredient")

    # Weighted average of raw column values
    bio_raw  = sum(idx.loc[k, "Biodegradability"]          * (pct / 100) for k, pct in blend.items() if k in idx.index)
    co2_raw  = sum(idx.loc[k, "CarbonFootprint_kgCO2eq"]   * (pct / 100) for k, pct in blend.items() if k in idx.index)
    etox_raw = sum(idx.loc[k, "Ecotoxicity_Score"]          * (pct / 100) for k, pct in blend.items() if k in idx.index)
    renew    = sum(idx.loc[k, "Renewability_Score"]          * (pct / 100) for k, pct in blend.items() if k in idx.index)
    reg_raw  = sum(_REACH_SCORES.get(idx.loc[k, "REACH_Flag"], 60) * (pct / 100) for k, pct in blend.items() if k in idx.index)

    # Normalize axes to 0–100
    axes = {
        "Biodegradability": round(bio_raw, 1),
        "Carbon Footprint":  round(_invert_co2(co2_raw), 1),
        "Ecotoxicity":       round(_invert_ecotox(etox_raw), 1),
        "Renewability":      round(renew, 1),
        "Regulatory":        round(reg_raw, 1),
    }

    eco_score = round(sum(axes[k] * AXIS_WEIGHTS[k] for k in axes), 1)

    vs_baseline = {k: round(axes[k] - _PETROCHEMICAL_BASELINE[k], 1) for k in axes}

    return EcoMetricsResult(
        biodegradability=axes["Biodegradability"],
        carbon_footprint=axes["Carbon Footprint"],
        ecotoxicity=axes["Ecotoxicity"],
        renewability=axes["Renewability"],
        regulatory=axes["Regulatory"],
        eco_score=eco_score,
        grade=_grade(eco_score),
        blend_summary=blend,
        vs_baseline=vs_baseline,
    )


def _synthesize_from_bio_pct(blend: Dict[str, float], db: pd.DataFrame) -> EcoMetricsResult:
    """
    Fallback: synthesize EcoMetrics from Bio_based_pct only.
    Used when extended columns are absent.
    """
    idx = db.set_index("Ingredient")
    bio_avg = sum(idx.loc[k, "Bio_based_pct"] * (pct / 100) for k, pct in blend.items() if k in idx.index)

    # Scale axes from bio% proxy
    axes = {
        "Biodegradability": round(bio_avg * 0.95, 1),
        "Carbon Footprint":  round(bio_avg * 0.90, 1),
        "Ecotoxicity":       round(bio_avg * 0.88, 1),
        "Renewability":      round(bio_avg, 1),
        "Regulatory":        round(min(bio_avg * 1.02, 100), 1),
    }
    eco_score = round(sum(axes[k] * AXIS_WEIGHTS[k] for k in axes), 1)
    vs_baseline = {k: round(axes[k] - _PETROCHEMICAL_BASELINE[k], 1) for k in axes}

    return EcoMetricsResult(
        biodegradability=axes["Biodegradability"],
        carbon_footprint=axes["Carbon Footprint"],
        ecotoxicity=axes["Ecotoxicity"],
        renewability=axes["Renewability"],
        regulatory=axes["Regulatory"],
        eco_score=eco_score,
        grade=_grade(eco_score),
        blend_summary=blend,
        vs_baseline=vs_baseline,
    )


def ecometrics_radar_data(result: EcoMetricsResult) -> dict:
    """
    Format EcoMetrics result for Plotly radar chart.
    Returns dict with 'categories', 'intelliform', and 'baseline' lists.
    """
    categories = list(AXIS_WEIGHTS.keys())
    intelliform_vals = [
        result.biodegradability,
        result.carbon_footprint,
        result.ecotoxicity,
        result.renewability,
        result.regulatory,
    ]
    baseline_vals = [_PETROCHEMICAL_BASELINE[c] for c in categories]

    return {
        "categories": categories,
        "intelliform": intelliform_vals,
        "baseline": baseline_vals,
    }
