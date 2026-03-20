"""
modules/carbon_credits.py
Carbon credit and circular economy calculator for IntelliForm v1.0.

Quantifies the carbon displacement value of switching from a
petrochemical baseline to an IntelliForm green formulation.

Methodology:
  - Petrochemical baseline: industry average INCI surfactant blend
    (Linear Alkylbenzene Sulfonate dominant, ~3.8 kgCO2eq/kg)
  - Green formulation: weighted average CarbonFootprint_kgCO2eq from DB
  - Carbon displacement = (baseline - green) * batch_kg
  - Credit value: Voluntary Carbon Market rate ($15-85/tonne CO2eq)
  - Circular economy score based on biodegradability + renewability

References:
  - GHG Protocol Product Standard
  - Voluntary Carbon Market Integrity Initiative (VCMI) 2023
  - ISO 14067:2018 Carbon footprint of products
"""
from dataclasses import dataclass
from typing import Dict
import pandas as pd


# Petrochemical baseline carbon footprint (kgCO2eq/kg)
# Based on: LAS surfactant industry average, BASF Product Carbon Footprint data
PETROCHEM_BASELINE_CO2 = 3.8

# Voluntary Carbon Market price range ($/tonne CO2eq) — 2026 VCM rates
VCM_PRICE_LOW  = 15.0   # floor — basic offsets
VCM_PRICE_MID  = 45.0   # mid — verified nature-based
VCM_PRICE_HIGH = 85.0   # premium — gold standard certified


@dataclass
class CarbonResult:
    # Per kg metrics
    green_co2_per_kg: float          # kgCO2eq/kg for this formulation
    baseline_co2_per_kg: float       # kgCO2eq/kg for petrochemical baseline
    displacement_per_kg: float       # kgCO2eq saved per kg produced

    # Per batch metrics (500 kg default)
    batch_kg: int
    co2_displaced_kg: float          # kg CO2eq displaced per batch
    co2_displaced_tonnes: float      # tonnes CO2eq displaced per batch

    # Credit value
    credit_value_low: float          # $ at VCM floor price
    credit_value_mid: float          # $ at VCM mid price
    credit_value_high: float         # $ at VCM premium price
    credits_per_batch: float         # number of carbon credits (1 credit = 1 tonne)

    # Annual projection (assuming 12 batches/year)
    annual_co2_tonnes: float
    annual_credit_value_mid: float

    # Circular economy
    circular_score: float            # 0-100
    circular_grade: str              # A+/A/B/C
    biodegradability_avg: float
    renewability_avg: float

    # Narrative
    summary: str


def calculate_carbon_credits(
    blend: Dict[str, float],
    db: pd.DataFrame,
    batch_kg: int = 500,
    batches_per_year: int = 12,
) -> CarbonResult:
    """
    Calculate carbon displacement and credit value for a formulation blend.
    """
    idx = db.set_index("Ingredient") if "Ingredient" in db.columns else db

    # ── Carbon footprint of green formulation ──
    green_co2 = 0.0
    bio_avg = 0.0
    renew_avg = 0.0
    covered_pct = 0.0

    for ing, pct in blend.items():
        if ing in idx.index:
            w = pct / 100
            try:
                co2 = float(idx.loc[ing, "CarbonFootprint_kgCO2eq"]) if "CarbonFootprint_kgCO2eq" in idx.columns else 1.5
                bio = float(idx.loc[ing, "Biodegradability"]) if "Biodegradability" in idx.columns else 85.0
                renew = float(idx.loc[ing, "Renewability_Score"]) if "Renewability_Score" in idx.columns else 80.0
                green_co2 += co2 * w
                bio_avg   += bio * w
                renew_avg += renew * w
                covered_pct += pct
            except Exception:
                green_co2 += 1.5 * w
                bio_avg   += 85.0 * w
                renew_avg += 80.0 * w

    # Fill any uncovered % with average
    if covered_pct < 99:
        remaining = (100 - covered_pct) / 100
        green_co2 += 1.5 * remaining
        bio_avg   += 85.0 * remaining
        renew_avg += 80.0 * remaining

    green_co2 = round(green_co2, 3)
    displacement_per_kg = round(PETROCHEM_BASELINE_CO2 - green_co2, 3)

    # ── Batch calculations ──
    co2_displaced_kg = round(displacement_per_kg * batch_kg, 1)
    co2_displaced_tonnes = round(co2_displaced_kg / 1000, 4)

    credit_value_low  = round(co2_displaced_tonnes * VCM_PRICE_LOW, 2)
    credit_value_mid  = round(co2_displaced_tonnes * VCM_PRICE_MID, 2)
    credit_value_high = round(co2_displaced_tonnes * VCM_PRICE_HIGH, 2)
    credits_per_batch = round(co2_displaced_tonnes, 4)

    # ── Annual projections ──
    annual_co2_tonnes = round(co2_displaced_tonnes * batches_per_year, 2)
    annual_credit_value_mid = round(annual_co2_tonnes * VCM_PRICE_MID, 2)

    # ── Circular economy score ──
    circular_score = round((bio_avg * 0.5) + (renew_avg * 0.5), 1)
    if circular_score >= 90: circular_grade = "A+"
    elif circular_score >= 80: circular_grade = "A"
    elif circular_score >= 70: circular_grade = "B"
    else: circular_grade = "C"

    # ── Summary narrative ──
    if displacement_per_kg > 0:
        summary = (
            f"This formulation displaces {displacement_per_kg:.2f} kgCO2eq per kg produced "
            f"vs petrochemical baseline. A {batch_kg}kg batch generates "
            f"{co2_displaced_tonnes:.3f} verified carbon credits worth "
            f"${credit_value_low:.0f}–${credit_value_high:.0f} on voluntary carbon markets. "
            f"At {batches_per_year} batches/year, annual carbon value is "
            f"~${annual_credit_value_mid:,.0f} (mid-market rate)."
        )
    else:
        summary = (
            f"This formulation has a carbon footprint of {green_co2:.2f} kgCO2eq/kg. "
            f"Consider substituting higher-carbon ingredients to generate positive carbon credits."
        )

    return CarbonResult(
        green_co2_per_kg=green_co2,
        baseline_co2_per_kg=PETROCHEM_BASELINE_CO2,
        displacement_per_kg=displacement_per_kg,
        batch_kg=batch_kg,
        co2_displaced_kg=co2_displaced_kg,
        co2_displaced_tonnes=co2_displaced_tonnes,
        credit_value_low=credit_value_low,
        credit_value_mid=credit_value_mid,
        credit_value_high=credit_value_high,
        credits_per_batch=credits_per_batch,
        annual_co2_tonnes=annual_co2_tonnes,
        annual_credit_value_mid=annual_credit_value_mid,
        circular_score=circular_score,
        circular_grade=circular_grade,
        biodegradability_avg=round(bio_avg, 1),
        renewability_avg=round(renew_avg, 1),
        summary=summary,
    )
