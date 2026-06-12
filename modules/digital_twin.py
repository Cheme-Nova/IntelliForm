"""
modules/digital_twin.py
Manufacturing Digital Twin — IntelliForm v1.6

A lightweight digital twin of the specialty-chemical manufacturing process.
Given a formulation blend and a target batch size, this module simulates:

  1. The process step timeline (charge, heat, mix, cool, QC hold) for the
     recommended manufacturing route
  2. Equipment sizing across lab → pilot → production scale tiers
  3. Energy consumption + CO2 (reuses carbon_passport energy intensities)
  4. Cycle time and daily throughput
  5. Scale-up risk factors (heat transfer, mixing time, foaming, etc.)
     and standard mitigations
  6. Yield-loss estimate vs ideal lab yield

This is NOT a CFD/process simulator — it is a transparent, rule-based twin
that gives formulators and toll manufacturers an immediate first-pass view
of "what does making this at scale actually look like" without requiring
plant data. As real ChemRich batch records accumulate, the heuristics below
can be replaced with fitted models without changing the public interface.

References:
  - Scale-up of mixing: Zlokarnik, "Scale-Up in Chemical Engineering" (2006)
  - Energy intensity by unit operation: see modules/carbon_passport.py
  - Digital twins for batch manufacturing: Rangaiah et al., 2020
"""
from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime
from typing import Dict, List, Optional
import math

import pandas as pd

from modules.carbon_passport import PROCESS_ENERGY, GRID_FACTORS
from modules.analytics import track


# ── Scale tiers ────────────────────────────────────────────────────────────────

SCALE_TIERS = {
    "lab":        {"label": "Lab",        "max_kg": 1,       "vessel": "Beaker / overhead stirrer (0.1–1 L)"},
    "bench":      {"label": "Bench",      "max_kg": 25,      "vessel": "Jacketed reactor (5–25 L)"},
    "pilot":      {"label": "Pilot",      "max_kg": 1000,    "vessel": "Pilot reactor / mixing tank (50 L–1 m³)"},
    "production": {"label": "Production", "max_kg": 50000,   "vessel": "Production reactor / tank farm (1–50 m³)"},
}


def _scale_tier(batch_kg: float) -> str:
    for tier, cfg in SCALE_TIERS.items():
        if batch_kg <= cfg["max_kg"]:
            return tier
    return "production"


# ── Process step templates by manufacturing process type ───────────────────────
# (duration_min is per-tonne-of-product baseline at "pilot" scale; scaled below)

PROCESS_STEP_TEMPLATES: Dict[str, List[Dict]] = {
    "blending": [
        {"name": "Charge raw materials",  "duration_min": 15, "temp_C": 20,  "ccp": False,
         "notes": "Verify each ingredient against COA before charging"},
        {"name": "Low-shear mixing",      "duration_min": 30, "temp_C": 25,  "ccp": True,
         "notes": "Mix until visually homogeneous; check for undissolved solids"},
        {"name": "QC sampling & hold",    "duration_min": 20, "temp_C": 25,  "ccp": True,
         "notes": "Pull sample for pH / viscosity / appearance before discharge"},
    ],
    "emulsification": [
        {"name": "Charge & heat aqueous phase",  "duration_min": 30, "temp_C": 70, "ccp": False,
         "notes": "Heat water phase to 70–75°C with slow agitation"},
        {"name": "Charge & heat oil phase",      "duration_min": 25, "temp_C": 75, "ccp": False,
         "notes": "Melt waxes/emulsifiers in separate vessel to 70–75°C"},
        {"name": "High-shear emulsification",    "duration_min": 20, "temp_C": 72, "ccp": True,
         "notes": "Add oil phase to aqueous phase under high shear (homogenizer)"},
        {"name": "Controlled cooling",           "duration_min": 60, "temp_C": 30, "ccp": True,
         "notes": "Cool with continued low-shear mixing to avoid phase separation"},
        {"name": "Add heat-sensitive actives",   "duration_min": 15, "temp_C": 35, "ccp": True,
         "notes": "Add preservatives, fragrance, actives below 40°C"},
        {"name": "QC sampling & hold",           "duration_min": 20, "temp_C": 25, "ccp": True,
         "notes": "Check viscosity, pH, droplet size / appearance"},
    ],
    "granulation": [
        {"name": "Dry blend powders",      "duration_min": 20, "temp_C": 22, "ccp": False,
         "notes": "Pre-blend API + excipients in high-shear granulator"},
        {"name": "Wet granulation",        "duration_min": 25, "temp_C": 35, "ccp": True,
         "notes": "Add binder solution gradually; monitor granule endpoint"},
        {"name": "Fluid-bed drying",       "duration_min": 90, "temp_C": 60, "ccp": True,
         "notes": "Dry to target LOD (typically <2%)"},
        {"name": "Milling / sizing",       "duration_min": 30, "temp_C": 22, "ccp": False,
         "notes": "Pass through oscillating granulator / sieve"},
        {"name": "Final blend & lubrication", "duration_min": 15, "temp_C": 22, "ccp": True,
         "notes": "Add lubricant, blend briefly to avoid over-lubrication"},
    ],
    "direct_compression": [
        {"name": "Dry blend (active + excipients)", "duration_min": 20, "temp_C": 22, "ccp": True,
         "notes": "Geometric dilution of low-dose actives; blend uniformity check"},
        {"name": "Lubrication blend",      "duration_min": 10, "temp_C": 22, "ccp": True,
         "notes": "Short blend time to avoid lubricant over-mixing"},
        {"name": "Compression",            "duration_min": 45, "temp_C": 22, "ccp": True,
         "notes": "Monitor tablet weight, hardness, thickness at intervals"},
    ],
    "spray_drying": [
        {"name": "Feed solution prep",     "duration_min": 40, "temp_C": 25, "ccp": False,
         "notes": "Dissolve/disperse solids in carrier solvent"},
        {"name": "Spray drying",           "duration_min": 60, "temp_C": 150, "ccp": True,
         "notes": "Inlet temp 120–180°C; monitor outlet temp and moisture"},
        {"name": "Powder collection & QC", "duration_min": 20, "temp_C": 25, "ccp": True,
         "notes": "Cyclone collection; check particle size distribution"},
    ],
    "extraction": [
        {"name": "Charge biomass + solvent", "duration_min": 30, "temp_C": 25, "ccp": False,
         "notes": "Load extraction vessel, verify solvent ratio"},
        {"name": "Extraction",              "duration_min": 120, "temp_C": 50, "ccp": True,
         "notes": "Maintain temperature and agitation per method validation"},
        {"name": "Filtration / separation", "duration_min": 45, "temp_C": 25, "ccp": True,
         "notes": "Remove spent biomass; clarify extract"},
        {"name": "Solvent recovery",        "duration_min": 60, "temp_C": 60, "ccp": False,
         "notes": "Distill solvent for recycle where possible"},
    ],
    "fermentation": [
        {"name": "Media prep & sterilization", "duration_min": 90, "temp_C": 121, "ccp": True,
         "notes": "Autoclave or in-situ sterilize fermentation media"},
        {"name": "Inoculation",             "duration_min": 15, "temp_C": 30, "ccp": True,
         "notes": "Aseptic transfer of starter culture"},
        {"name": "Fermentation",            "duration_min": 1440, "temp_C": 30, "ccp": True,
         "notes": "Monitor pH, DO, temperature continuously (~24h cycle)"},
        {"name": "Harvest & downstream",    "duration_min": 120, "temp_C": 25, "ccp": True,
         "notes": "Centrifuge / filter biomass, recover product stream"},
    ],
    "default": [
        {"name": "Charge raw materials",  "duration_min": 20, "temp_C": 22, "ccp": False,
         "notes": "Verify materials against formulation sheet"},
        {"name": "Process / mix",         "duration_min": 45, "temp_C": 30, "ccp": True,
         "notes": "Process per standard operating procedure for this product class"},
        {"name": "QC sampling & hold",    "duration_min": 20, "temp_C": 25, "ccp": True,
         "notes": "Release testing before discharge/packaging"},
    ],
}

# Map free-form manufacturing_route strings (e.g. from optimizer/pharma) to a
# PROCESS_STEP_TEMPLATES / PROCESS_ENERGY key
ROUTE_TO_PROCESS = {
    "direct compression": "direct_compression",
    "wet granulation": "granulation",
    "dry granulation": "granulation",
    "roller compaction": "granulation",
    "matrix tablet": "granulation",
    "cream": "emulsification",
    "lotion": "emulsification",
    "emulsion": "emulsification",
    "aseptic": "spray_drying",
}


def _resolve_process(manufacturing_process: Optional[str],
                      manufacturing_route: Optional[str]) -> str:
    """Resolve a process-step-template key from explicit process or route hints."""
    if manufacturing_process and manufacturing_process in PROCESS_STEP_TEMPLATES:
        return manufacturing_process
    if manufacturing_route:
        route_l = manufacturing_route.lower()
        for key, proc in ROUTE_TO_PROCESS.items():
            if key in route_l:
                return proc
    return "blending"


# ── Scale-up risk factors ────────────────────────────────────────────────────

SCALE_UP_RISKS: Dict[str, Dict] = {
    "mixing_time": {
        "factor": "Mixing time scale-up",
        "trigger_processes": {"blending", "emulsification", "granulation"},
        "mitigation": "Mixing time generally increases with vessel diameter — "
                       "verify homogeneity with in-process sampling at 3 vessel heights, "
                       "not just impeller-tip speed matching.",
    },
    "heat_transfer": {
        "factor": "Heat transfer rate",
        "trigger_processes": {"emulsification", "granulation", "spray_drying", "fermentation"},
        "mitigation": "Surface-area-to-volume ratio drops sharply at scale — heating/cooling "
                       "steps will take longer than lab timing suggests. Add temperature "
                       "ramp profiling to the batch record at pilot scale before production.",
    },
    "foaming": {
        "factor": "Foaming / air entrainment",
        "trigger_processes": {"emulsification"},
        "mitigation": "High-shear steps entrain more air at larger vessel scale. "
                       "Consider vacuum deaeration step or anti-foam at pilot trial.",
    },
    "granule_attrition": {
        "factor": "Granule attrition / particle size shift",
        "trigger_processes": {"granulation"},
        "mitigation": "Fluid-bed dryer airflow and impeller tip speed both change "
                       "particle size distribution at scale — run a bridging study "
                       "comparing pilot and production PSD before validation.",
    },
    "drying_uniformity": {
        "factor": "Drying / moisture uniformity",
        "trigger_processes": {"spray_drying", "granulation"},
        "mitigation": "Larger dryers create non-uniform residence time — sample "
                       "multiple zones for moisture content (LOD) during scale-up runs.",
    },
    "sterility_assurance": {
        "factor": "Sterility / contamination control",
        "trigger_processes": {"fermentation", "spray_drying"},
        "mitigation": "Aseptic processing risk increases with batch duration and "
                       "open transfers — formal contamination risk assessment required "
                       "before production-scale runs.",
    },
}


# ── Result schemas ───────────────────────────────────────────────────────────

@dataclass
class ProcessStep:
    order: int
    name: str
    duration_min: float
    temperature_C: float
    critical_control_point: bool
    notes: str


@dataclass
class ScaleUpRisk:
    factor: str
    risk_level: str       # Low / Medium / High
    mitigation: str


@dataclass
class DigitalTwinResult:
    batch_kg: float
    scale_tier: str
    scale_label: str
    recommended_vessel: str
    manufacturing_process: str
    process_steps: List[ProcessStep]
    total_cycle_time_min: float
    total_cycle_time_hr: float
    batches_per_day: float
    energy_kwh_per_batch: float
    energy_cost_usd_per_batch: float
    co2_kg_per_batch: float
    co2_kg_per_kg: float
    yield_pct: float
    output_kg: float
    scale_up_risks: List[ScaleUpRisk]
    notes: str


@dataclass
class ScaleComparison:
    tiers: List[DigitalTwinResult]
    summary: str


# ── Core simulation ───────────────────────────────────────────────────────────

# Approximate industrial electricity price by region (USD/kWh) — for cost estimates only
ENERGY_PRICE_USD_PER_KWH = {
    "US": 0.08, "EU": 0.18, "UK": 0.20, "China": 0.09, "India": 0.08, "Global": 0.12,
}

# Baseline yield by scale tier — small batches lose proportionally more to
# transfer/hold-up losses; larger batches benefit from fixed-loss dilution.
BASE_YIELD_BY_TIER = {
    "lab": 0.92, "bench": 0.95, "pilot": 0.97, "production": 0.985,
}


def _mixing_time_scale_factor(batch_kg: float, reference_kg: float = 50.0) -> float:
    """
    Heuristic mixing-time scale-up: time scales roughly with the cube root of
    the volume ratio (geometric similarity, constant tip speed assumption),
    bounded to avoid runaway values for very large/small batches.
    """
    ratio = max(batch_kg, 0.01) / reference_kg
    factor = ratio ** (1.0 / 3.0)
    return float(min(max(factor, 0.5), 4.0))


def simulate_manufacturing(
    blend: Dict[str, float],
    db: pd.DataFrame,
    batch_kg: float = 500.0,
    manufacturing_process: Optional[str] = None,
    manufacturing_route: Optional[str] = None,
    vertical: str = "all",
    grid_region: str = "US",
) -> DigitalTwinResult:
    """
    Simulate the manufacturing process for a formulation at a given batch size.

    Args:
        blend: {ingredient: percentage} — used for risk heuristics (oil phase load etc.)
        db: ingredient database (used for function lookups)
        batch_kg: target batch size in kg
        manufacturing_process: explicit process key (e.g. "emulsification");
            if not provided, resolved from manufacturing_route or vertical
        manufacturing_route: free-form route string (from optimizer/pharma modules)
        vertical: industry vertical — used as a fallback hint for process selection
        grid_region: electricity grid region for CO2/cost (see carbon_passport.GRID_FACTORS)

    Returns:
        DigitalTwinResult with process timeline, energy/CO2, yield, and scale-up risks.
    """
    process = _resolve_process(manufacturing_process, manufacturing_route)
    if process == "blending" and not manufacturing_process and not manufacturing_route:
        # Vertical-based fallback when nothing else is known
        vertical_defaults = {
            "personal_care": "emulsification",
            "pharmaceutical": "direct_compression",
            "food_beverage": "blending",
            "agriculture": "blending",
            "industrial": "blending",
            "home_care": "blending",
            "cosmetics": "emulsification",
        }
        process = vertical_defaults.get(vertical, "blending")

    tier = _scale_tier(batch_kg)
    tier_cfg = SCALE_TIERS[tier]

    templates = PROCESS_STEP_TEMPLATES.get(process, PROCESS_STEP_TEMPLATES["default"])
    time_factor = _mixing_time_scale_factor(batch_kg)

    steps: List[ProcessStep] = []
    total_min = 0.0
    for i, t in enumerate(templates, start=1):
        # Only mixing/heating/cooling-type CCP steps scale with batch size;
        # fixed steps (charging, sampling) scale much less
        scale = time_factor if t["ccp"] else min(time_factor, 1.5)
        duration = round(t["duration_min"] * scale, 1)
        steps.append(ProcessStep(
            order=i, name=t["name"], duration_min=duration,
            temperature_C=t["temp_C"], critical_control_point=t["ccp"],
            notes=t["notes"],
        ))
        total_min += duration

    total_hr = total_min / 60.0
    # Assume one shift (8h effective) per batch slot; allow for changeover
    batches_per_day = round(max(24.0 / (total_hr + 1.0), 0.1), 2)

    # ── Energy & CO2 (reuse carbon_passport intensities) ────────────────────
    energy_kwh_per_kg = PROCESS_ENERGY.get(process, PROCESS_ENERGY["default"])
    energy_kwh = energy_kwh_per_kg * batch_kg
    grid_factor = GRID_FACTORS.get(grid_region, GRID_FACTORS["Global"])
    co2_kg = energy_kwh * grid_factor
    price = ENERGY_PRICE_USD_PER_KWH.get(grid_region, ENERGY_PRICE_USD_PER_KWH["Global"])
    energy_cost = energy_kwh * price

    # ── Yield ────────────────────────────────────────────────────────────────
    yield_pct = BASE_YIELD_BY_TIER.get(tier, 0.95) * 100
    # Fermentation / extraction carry inherent extra losses
    if process in ("fermentation", "extraction"):
        yield_pct -= 5
    output_kg = round(batch_kg * yield_pct / 100, 2)

    # ── Scale-up risks ───────────────────────────────────────────────────────
    risks: List[ScaleUpRisk] = []
    for risk in SCALE_UP_RISKS.values():
        if process not in risk["trigger_processes"]:
            continue
        if tier == "lab":
            continue  # not yet scaling
        level = "High" if tier == "production" else "Medium"
        risks.append(ScaleUpRisk(
            factor=risk["factor"], risk_level=level, mitigation=risk["mitigation"],
        ))

    # Oil-phase-specific note for emulsification at scale
    oil_pct = sum(v for ing, v in blend.items()
                   if any(o in ing.lower() for o in ["oil", "ester", "wax", "butter"]))
    if process == "emulsification" and oil_pct > 30 and tier != "lab":
        risks.append(ScaleUpRisk(
            factor="High oil-phase load at scale",
            risk_level="Medium",
            mitigation=f"Oil phase is {oil_pct:.0f}% of formula — confirm homogenizer "
                       f"shear rate is sufficient at production vessel diameter; "
                       f"consider in-line homogenization for batches >1,000 kg.",
        ))

    notes = (
        f"Resolved manufacturing process: {process.replace('_', ' ')}. "
        f"Scale tier: {tier_cfg['label']} ({tier_cfg['vessel']}). "
        f"Cycle time scales with batch size via a cube-root mixing-time heuristic "
        f"(reference batch: 50 kg)."
    )

    track("digital_twin_simulated", {
        "batch_kg": batch_kg, "process": process, "tier": tier, "vertical": vertical,
    })

    return DigitalTwinResult(
        batch_kg=batch_kg,
        scale_tier=tier,
        scale_label=tier_cfg["label"],
        recommended_vessel=tier_cfg["vessel"],
        manufacturing_process=process,
        process_steps=steps,
        total_cycle_time_min=round(total_min, 1),
        total_cycle_time_hr=round(total_hr, 2),
        batches_per_day=batches_per_day,
        energy_kwh_per_batch=round(energy_kwh, 2),
        energy_cost_usd_per_batch=round(energy_cost, 2),
        co2_kg_per_batch=round(co2_kg, 2),
        co2_kg_per_kg=round(co2_kg / batch_kg, 4) if batch_kg else 0.0,
        yield_pct=round(yield_pct, 1),
        output_kg=output_kg,
        scale_up_risks=risks,
        notes=notes,
    )


def compare_scales(
    blend: Dict[str, float],
    db: pd.DataFrame,
    manufacturing_process: Optional[str] = None,
    manufacturing_route: Optional[str] = None,
    vertical: str = "all",
    grid_region: str = "US",
    scales_kg: Optional[List[float]] = None,
) -> ScaleComparison:
    """
    Run the digital twin across lab → pilot → production scale points so users
    can see how cycle time, energy, and risk profile change with batch size.
    """
    if scales_kg is None:
        scales_kg = [1, 50, 1000, 10000]

    results = [
        simulate_manufacturing(
            blend, db, batch_kg=kg,
            manufacturing_process=manufacturing_process,
            manufacturing_route=manufacturing_route,
            vertical=vertical, grid_region=grid_region,
        )
        for kg in scales_kg
    ]

    lab, prod = results[0], results[-1]
    energy_ratio = (prod.co2_kg_per_kg / lab.co2_kg_per_kg) if lab.co2_kg_per_kg else 1.0
    summary = (
        f"From {lab.batch_kg:g} kg ({lab.scale_label}) to {prod.batch_kg:g} kg "
        f"({prod.scale_label}): cycle time {lab.total_cycle_time_hr:.1f}h → "
        f"{prod.total_cycle_time_hr:.1f}h, CO2/kg {'increases' if energy_ratio > 1 else 'decreases'} "
        f"by {abs(energy_ratio - 1) * 100:.0f}%, yield {lab.yield_pct:.1f}% → {prod.yield_pct:.1f}%. "
        f"{len(prod.scale_up_risks)} scale-up risk factor(s) identified at production scale."
    )

    return ScaleComparison(tiers=results, summary=summary)
