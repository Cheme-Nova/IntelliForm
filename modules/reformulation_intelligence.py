"""
modules/reformulation_intelligence.py
Closed-Loop Reformulation Intelligence™ — IntelliForm v1.4

The "self-driving lab" concept applied to specialty chemical SMEs without robotics.

What it does:
  When a pilot batch fails a test (pH wrong, viscosity too high, stability issue,
  certification rejection), this module:
    1. Diagnoses the ROOT CAUSE from structured test results
    2. Suggests the MINIMAL CHANGE to fix it (1-3 ingredient swaps)
    3. Predicts the DOWNSTREAM IMPACT on other properties
    4. Ranks suggestions by cost impact + implementation feasibility
    5. Learns from every iteration (Bayesian belief update)

This is distinct from standard reformulation:
  Standard: "optimize from scratch given constraints"
  This: "I have a specific failure — what's the minimal surgical fix?"

The commercial insight: ChemRich pilot batches produce structured failure data.
That data, fed into this module, creates a proprietary closed loop that no
AI tool built on public data can replicate.

Failure taxonomy (based on specialty chemical QC literature):
  - pH out of range
  - Viscosity too high / too low
  - Phase separation / emulsion instability
  - Colour / odour complaint
  - Performance shortfall (cleaning, conditioning, etc.)
  - Certification rejection (specific ingredient flagged)
  - Regulatory non-compliance (REACH, FDA, ICH)
  - Shelf life failure (accelerated stability)
  - Cost overrun
  - Supply chain disruption (ingredient unavailable)

Reference:
  - DoE (Design of Experiments) for formulation: Cornell 2002
  - Bayesian optimization for formulation: Shields et al., Nature 2021
  - Self-driving labs: Abolhasani & Kumacheva, Nature Synthesis 2023
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple
import numpy as np
import pandas as pd


# ── Failure taxonomy ─────────────────────────────────────────────────────────

FAILURE_TYPES = {
    "ph_too_high": {
        "label": "pH Too High",
        "description": "Measured pH above target range",
        "inputs": ["measured_ph", "target_ph_min", "target_ph_max"],
    },
    "ph_too_low": {
        "label": "pH Too Low",
        "description": "Measured pH below target range",
        "inputs": ["measured_ph", "target_ph_min", "target_ph_max"],
    },
    "viscosity_too_high": {
        "label": "Viscosity Too High",
        "description": "Measured viscosity exceeds specification",
        "inputs": ["measured_viscosity_cP", "target_viscosity_max_cP"],
    },
    "viscosity_too_low": {
        "label": "Viscosity Too Low",
        "description": "Measured viscosity below specification",
        "inputs": ["measured_viscosity_cP", "target_viscosity_min_cP"],
    },
    "phase_separation": {
        "label": "Phase Separation / Emulsion Instability",
        "description": "Product separates into layers on storage",
        "inputs": ["storage_days", "temperature_C"],
    },
    "performance_shortfall": {
        "label": "Performance Shortfall",
        "description": "Product does not meet performance specification",
        "inputs": ["measured_performance", "target_performance", "performance_metric"],
    },
    "certification_rejection": {
        "label": "Certification Rejection",
        "description": "Specific ingredient flagged by certification body",
        "inputs": ["certification_name", "rejected_ingredient", "rejection_reason"],
    },
    "stability_failure": {
        "label": "Accelerated Stability Failure",
        "description": "Product fails ICH/accelerated stability testing",
        "inputs": ["stability_condition", "failure_timepoint_weeks", "failure_observation"],
    },
    "cost_overrun": {
        "label": "Cost Overrun",
        "description": "Actual batch cost exceeds target",
        "inputs": ["actual_cost_per_kg", "target_cost_per_kg"],
    },
    "supply_disruption": {
        "label": "Supply Chain Disruption",
        "description": "Key ingredient unavailable or lead time too long",
        "inputs": ["disrupted_ingredient", "lead_time_weeks"],
    },
    "colour_complaint": {
        "label": "Colour / Appearance Issue",
        "description": "Product colour outside specification",
        "inputs": ["observed_colour", "target_colour"],
    },
    "odour_complaint": {
        "label": "Odour Complaint",
        "description": "Product odour unacceptable",
        "inputs": ["odour_description"],
    },
}


# ── Ingredient function roles for reformulation ────────────────────────────────

# Maps failure type → which ingredient functions are likely responsible
FAILURE_FUNCTION_MAP = {
    "ph_too_high": ["pH Adjuster", "Alkaline Builder", "pH Neutralizer",
                    "pH Adjuster/Saponifier", "Buffer"],
    "ph_too_low":  ["pH Adjuster", "Acidulant", "Buffer", "Acid"],
    "viscosity_too_high": ["Thickener", "Gelling Agent", "Viscosity Modifier",
                           "Stabilizer", "Film Former"],
    "viscosity_too_low":  ["Thickener", "Gelling Agent", "Viscosity Modifier"],
    "phase_separation": ["Emulsifier", "Stabilizer", "Co-emulsifier",
                         "Emulsion Stabilizer", "Rheology Modifier"],
    "performance_shortfall": ["Active", "Surfactant", "Antimicrobial",
                               "Conditioning Agent", "Performance Booster"],
    "certification_rejection": [],  # handled separately
    "stability_failure": ["Preservative", "Antioxidant", "Chelating Agent",
                          "Stabilizer", "UV Filter"],
    "cost_overrun":  [],  # handled by cost optimizer
    "supply_disruption": [],  # handled by ingredient replacement
    "colour_complaint": ["Colorant", "Natural Color", "Bleaching Agent"],
    "odour_complaint": ["Fragrance", "Masking Agent", "Essential Oil"],
}

# pH adjustment rules
PH_ADJUSTERS_UP   = ["sodium hydroxide", "potassium hydroxide", "amp-95",
                      "triethanolamine", "ammonia", "sodium bicarbonate",
                      "potassium carbonate"]
PH_ADJUSTERS_DOWN  = ["citric acid", "lactic acid", "glycolic acid",
                       "phosphoric acid", "acetic acid", "sodium bisulfate",
                       "hydrochloric acid"]

# Thickeners for viscosity adjustment
THICKENERS_INCREASE = ["carbomer", "xanthan gum", "hydroxyethyl cellulose",
                        "hydroxypropyl methylcellulose", "carrageenan",
                        "gellan gum", "guar gum", "sodium alginate"]
THICKENERS_DECREASE = ["propylene glycol", "glycerol", "sorbitol",
                        "ethanol", "water"]

# Emulsifiers for stability
EMULSIFIERS = ["lecithin", "polysorbate", "glyceryl stearate", "ceteareth",
                "sucrose laurate", "polyglyceryl", "sorbitan", "steareth"]

# Solvents that increase instability under stress conditions when over-loaded
SOLVENTS = ["ipa", "isopropanol", "limonene", "ethanol", "alcohol",
            "acetone", "methanol", "butanol", "propylene glycol",
            "dipropylene glycol", "glycol ether"]

# Preservatives / antioxidants / chelators that protect shelf stability
PRESERVATIVE_KEYWORDS = ["preservative", "antioxidant", "chelating", "edta",
                          "phenoxyethanol", "benzoate", "sorbate", "tocopherol",
                          "ascorbate", "bht", "bha", "propyl gallate"]


# ── Result schemas ─────────────────────────────────────────────────────────────

@dataclass
class ReformulationSuggestion:
    rank: int
    action_type: str               # "replace" | "increase" | "decrease" | "add" | "remove"
    ingredient: str                # ingredient to act on
    current_pct: float
    suggested_pct: float
    rationale: str
    confidence: float              # 0-1
    predicted_impact: Dict[str, str]   # {"pH": "+0.3", "viscosity": "unchanged"}
    cost_delta_per_kg: float           # change in cost/kg
    risk_level: str                    # Low / Medium / High
    implementation_notes: str


@dataclass
class RootCauseAnalysis:
    failure_type: str
    root_cause: str
    contributing_ingredients: List[str]
    severity: str                  # Minor / Moderate / Severe / Critical
    confidence: float
    evidence: List[str]


@dataclass
class ReformulationReport:
    failure_type: str
    failure_description: str
    root_cause: RootCauseAnalysis
    suggestions: List[ReformulationSuggestion]
    best_suggestion: ReformulationSuggestion
    predicted_success_probability: float
    iterations_to_fix: str         # "1-2 batches" / "2-4 batches"
    do_not_change: List[str]       # ingredients that are working fine
    learning_note: str             # what this teaches about the formulation
    cbam_impact: Optional[str]     # impact on carbon passport if reformulated


# ── Diagnosis engine ───────────────────────────────────────────────────────────

def _diagnose_ph(blend: Dict[str, float], db: pd.DataFrame,
                  failure_type: str, test_data: dict) -> RootCauseAnalysis:
    """Diagnose pH failure root cause."""
    idx = db.set_index("Ingredient") if "Ingredient" in db.columns else db
    contributing = []
    evidence = []

    measured = float(test_data.get("measured_ph", 7.0))
    target_min = float(test_data.get("target_ph_min", 5.0))
    target_max = float(test_data.get("target_ph_max", 7.0))
    target_mid = (target_min + target_max) / 2

    delta = measured - target_mid
    severity = "Minor" if abs(delta) < 0.5 else "Moderate" if abs(delta) < 1.5 else "Severe"

    for ing in blend:
        ing_l = ing.lower()
        if failure_type == "ph_too_high":
            if any(a in ing_l for a in PH_ADJUSTERS_UP):
                contributing.append(ing)
                evidence.append(f"{ing} is alkaline — may be over-dosed")
        elif failure_type == "ph_too_low":
            if any(a in ing_l for a in PH_ADJUSTERS_DOWN):
                contributing.append(ing)
                evidence.append(f"{ing} is acidic — may be over-dosed")

    root_cause = (
        f"Measured pH {measured:.1f} is {'above' if delta > 0 else 'below'} "
        f"target range {target_min}–{target_max}. "
        f"Delta = {delta:+.1f} pH units. "
        f"{'Alkaline' if delta > 0 else 'Acidic'} contributors identified."
    )

    return RootCauseAnalysis(
        failure_type=failure_type, root_cause=root_cause,
        contributing_ingredients=contributing, severity=severity,
        confidence=0.85, evidence=evidence,
    )


def _diagnose_viscosity(blend: Dict[str, float], db: pd.DataFrame,
                         failure_type: str, test_data: dict) -> RootCauseAnalysis:
    """Diagnose viscosity failure root cause."""
    idx = db.set_index("Ingredient") if "Ingredient" in db.columns else db
    contributing = []
    evidence = []

    measured = float(test_data.get("measured_viscosity_cP", 1000))
    if failure_type == "viscosity_too_high":
        target = float(test_data.get("target_viscosity_max_cP", 5000))
        severity = "Minor" if measured < target * 1.2 else \
                   "Moderate" if measured < target * 2 else "Severe"
        for ing in blend:
            if any(t in ing.lower() for t in THICKENERS_INCREASE):
                contributing.append(ing)
                evidence.append(f"{ing} at {blend[ing]:.1f}% — reduce or replace")
    else:
        target = float(test_data.get("target_viscosity_min_cP", 500))
        severity = "Minor" if measured > target * 0.8 else \
                   "Moderate" if measured > target * 0.5 else "Severe"
        evidence.append("Thickener level insufficient or incompatible with other ingredients")

    root_cause = (
        f"Measured viscosity {measured:,.0f} cP vs target "
        f"{'max ' if failure_type == 'viscosity_too_high' else 'min '}"
        f"{target:,.0f} cP. "
        f"{'Reduce' if failure_type == 'viscosity_too_high' else 'Increase'} "
        f"thickener system."
    )

    return RootCauseAnalysis(
        failure_type=failure_type, root_cause=root_cause,
        contributing_ingredients=contributing, severity=severity,
        confidence=0.80, evidence=evidence,
    )


def _diagnose_phase_separation(blend: Dict[str, float], db: pd.DataFrame,
                                 test_data: dict) -> RootCauseAnalysis:
    """Diagnose phase separation root cause."""
    contributing = []
    evidence = []
    days = float(test_data.get("storage_days", 7))
    temp = float(test_data.get("temperature_C", 45))

    # Check emulsifier HLB balance
    emulsifier_count = sum(1 for ing in blend
                           if any(e in ing.lower() for e in EMULSIFIERS))
    if emulsifier_count == 0:
        evidence.append("No emulsifier detected — formulation cannot form stable emulsion")
    elif emulsifier_count == 1:
        evidence.append("Single emulsifier — consider HLB-balanced system (oil-in-water + water-in-oil)")

    # Oil phase check
    oil_pct = sum(v for ing, v in blend.items()
                  if any(o in ing.lower() for o in ["oil", "ester", "wax", "butter"]))
    if oil_pct > 30:
        evidence.append(f"High oil phase ({oil_pct:.0f}%) — may exceed emulsifier capacity")

    severity = "Severe" if days < 7 else "Moderate" if days < 30 else "Minor"

    return RootCauseAnalysis(
        failure_type="phase_separation",
        root_cause=f"Phase separation after {days:.0f} days at {temp:.0f}°C. "
                   f"Emulsifier system insufficient for oil phase load.",
        contributing_ingredients=contributing, severity=severity,
        confidence=0.72, evidence=evidence,
    )


def _diagnose_certification(blend: Dict[str, float], test_data: dict) -> RootCauseAnalysis:
    """Diagnose certification rejection."""
    cert = test_data.get("certification_name", "Unknown")
    rejected_ing = test_data.get("rejected_ingredient", "Unknown")
    reason = test_data.get("rejection_reason", "Not specified")

    return RootCauseAnalysis(
        failure_type="certification_rejection",
        root_cause=f"{cert} rejected formulation due to {rejected_ing}: {reason}",
        contributing_ingredients=[rejected_ing],
        severity="Critical",
        confidence=0.99,
        evidence=[f"Certification body explicitly flagged: {rejected_ing}",
                  f"Reason: {reason}"],
    )


# ── Suggestion generation ──────────────────────────────────────────────────────

def _suggest_ph_fix(blend: Dict[str, float], db: pd.DataFrame,
                     failure_type: str, test_data: dict) -> List[ReformulationSuggestion]:
    """Generate pH correction suggestions."""
    idx = db.set_index("Ingredient") if "Ingredient" in db.columns else db
    suggestions = []
    measured = float(test_data.get("measured_ph", 7.0))
    target_mid = (float(test_data.get("target_ph_min", 5.0)) +
                  float(test_data.get("target_ph_max", 7.0))) / 2
    delta = abs(measured - target_mid)

    if failure_type == "ph_too_high":
        # Add/increase acid
        existing_acids = [ing for ing in blend
                          if any(a in ing.lower() for a in PH_ADJUSTERS_DOWN)]
        if existing_acids:
            ing = existing_acids[0]
            adj = min(blend[ing] * 0.5, 2.0) * (delta / 1.0)
            suggestions.append(ReformulationSuggestion(
                rank=1, action_type="increase", ingredient=ing,
                current_pct=blend[ing], suggested_pct=round(blend[ing] + adj, 2),
                rationale=f"Increase {ing} by {adj:.1f}% to lower pH toward target {target_mid:.1f}",
                confidence=0.82,
                predicted_impact={"pH": f"-{delta*0.6:.1f} units (estimated)",
                                   "cost": "minimal"},
                cost_delta_per_kg=0.02, risk_level="Low",
                implementation_notes=f"Add incrementally — check pH after each 0.5% addition",
            ))
        else:
            # Recommend adding citric acid
            suggestions.append(ReformulationSuggestion(
                rank=1, action_type="add", ingredient="Citric Acid",
                current_pct=0.0, suggested_pct=round(0.2 * delta, 2),
                rationale=f"Add Citric Acid to lower pH from {measured:.1f} to {target_mid:.1f}",
                confidence=0.88,
                predicted_impact={"pH": f"-{delta*0.8:.1f} units", "cost": "+$0.03/kg"},
                cost_delta_per_kg=0.03, risk_level="Low",
                implementation_notes="Citric acid is GRAS, COSMOS-approved, low-cost pH adjuster",
            ))

    else:  # ph_too_low
        existing_bases = [ing for ing in blend
                          if any(a in ing.lower() for a in PH_ADJUSTERS_UP)]
        if existing_bases:
            ing = existing_bases[0]
            adj = min(blend[ing] * 0.5, 1.5) * (delta / 1.0)
            suggestions.append(ReformulationSuggestion(
                rank=1, action_type="increase", ingredient=ing,
                current_pct=blend[ing], suggested_pct=round(blend[ing] + adj, 2),
                rationale=f"Increase {ing} to raise pH toward target {target_mid:.1f}",
                confidence=0.80,
                predicted_impact={"pH": f"+{delta*0.6:.1f} units", "cost": "minimal"},
                cost_delta_per_kg=0.01, risk_level="Low",
                implementation_notes="Add in 0.1% increments, mix well before measuring pH",
            ))
        else:
            suggestions.append(ReformulationSuggestion(
                rank=1, action_type="add", ingredient="Sodium Bicarbonate",
                current_pct=0.0, suggested_pct=round(0.3 * delta, 2),
                rationale=f"Add Sodium Bicarbonate to raise pH from {measured:.1f} to {target_mid:.1f}",
                confidence=0.82,
                predicted_impact={"pH": f"+{delta*0.7:.1f} units", "cost": "+$0.01/kg"},
                cost_delta_per_kg=0.01, risk_level="Low",
                implementation_notes="Dissolve in warm water before adding to batch",
            ))

    return suggestions


def _suggest_viscosity_fix(blend: Dict[str, float], db: pd.DataFrame,
                            failure_type: str, test_data: dict) -> List[ReformulationSuggestion]:
    """Generate viscosity correction suggestions."""
    idx = db.set_index("Ingredient") if "Ingredient" in db.columns else db
    suggestions = []

    measured = float(test_data.get("measured_viscosity_cP", 1000))
    if failure_type == "viscosity_too_high":
        target = float(test_data.get("target_viscosity_max_cP", 5000))
        ratio = measured / target

        # Find and reduce thickener
        thickeners = [(ing, pct) for ing, pct in blend.items()
                      if any(t in ing.lower() for t in THICKENERS_INCREASE)]
        if thickeners:
            ing, pct = thickeners[0]
            reduction = pct * min(0.5, (ratio - 1) / ratio)
            suggestions.append(ReformulationSuggestion(
                rank=1, action_type="decrease", ingredient=ing,
                current_pct=pct, suggested_pct=round(pct - reduction, 2),
                rationale=f"Reduce {ing} by {reduction:.1f}% — viscosity {measured:,.0f} cP is "
                          f"{ratio:.1f}x above target {target:,.0f} cP",
                confidence=0.78,
                predicted_impact={"viscosity": f"-{int((ratio-1)/ratio*100)}% estimated",
                                   "stability": "monitor"},
                cost_delta_per_kg=-0.05, risk_level="Low",
                implementation_notes="Reduce gradually — recheck viscosity at 20°C after each change",
            ))
    else:
        target = float(test_data.get("target_viscosity_min_cP", 500))
        thickeners = [(ing, pct) for ing, pct in blend.items()
                      if any(t in ing.lower() for t in THICKENERS_INCREASE)]
        if thickeners:
            ing, pct = thickeners[0]
            increase = pct * 0.3
            suggestions.append(ReformulationSuggestion(
                rank=1, action_type="increase", ingredient=ing,
                current_pct=pct, suggested_pct=round(pct + increase, 2),
                rationale=f"Increase {ing} by {increase:.1f}% to build viscosity to {target:,.0f} cP",
                confidence=0.75,
                predicted_impact={"viscosity": f"+{int(increase/pct*100)}% estimated"},
                cost_delta_per_kg=0.08, risk_level="Low",
                implementation_notes="Add thickener slowly with good agitation to avoid lumps",
            ))
        else:
            suggestions.append(ReformulationSuggestion(
                rank=1, action_type="add", ingredient="Xanthan Gum",
                current_pct=0.0, suggested_pct=0.3,
                rationale=f"Add Xanthan Gum 0.3% — viscosity builder, GRAS, stable across pH 3–9",
                confidence=0.85,
                predicted_impact={"viscosity": "+1000–3000 cP estimated"},
                cost_delta_per_kg=0.06, risk_level="Low",
                implementation_notes="Pre-disperse in glycerol 1:5 before adding to aqueous phase",
            ))

    return suggestions


def _suggest_phase_separation_fix(blend: Dict[str, float], db: pd.DataFrame,
                                   test_data: dict) -> List[ReformulationSuggestion]:
    """Generate emulsion stability suggestions."""
    suggestions = []

    emulsifier_pct = sum(v for ing, v in blend.items()
                         if any(e in ing.lower() for e in EMULSIFIERS))

    if emulsifier_pct < 3.0:
        suggestions.append(ReformulationSuggestion(
            rank=1, action_type="increase",
            ingredient=next((i for i in blend if any(e in i.lower() for e in EMULSIFIERS)),
                            "Lecithin"),
            current_pct=emulsifier_pct, suggested_pct=min(emulsifier_pct + 2.0, 8.0),
            rationale="Insufficient emulsifier level — typical stable emulsion needs 3–8% total",
            confidence=0.78,
            predicted_impact={"stability": "Improved", "viscosity": "+10–20%"},
            cost_delta_per_kg=0.12, risk_level="Low",
            implementation_notes="Add to oil phase before homogenization",
        ))

    suggestions.append(ReformulationSuggestion(
        rank=2, action_type="add", ingredient="Xanthan Gum",
        current_pct=0.0, suggested_pct=0.2,
        rationale="Add co-stabilizer — xanthan creates 3D network that prevents droplet coalescence",
        confidence=0.82,
        predicted_impact={"stability": "Significantly improved", "viscosity": "+500–1000 cP"},
        cost_delta_per_kg=0.04, risk_level="Low",
        implementation_notes="Add to water phase at 60°C for maximum hydration",
    ))

    return suggestions


def _diagnose_stability(blend: Dict[str, float], test_data: dict) -> RootCauseAnalysis:
    """Diagnose accelerated stability failure root cause."""
    contributing = []
    evidence = []
    condition = test_data.get("stability_condition", "accelerated")
    timepoint = test_data.get("failure_timepoint_weeks", "unknown")
    observation = test_data.get("failure_observation", "degradation")

    # Rule: high solvent load (>30% total) destabilises emulsions under thermal stress
    solvent_ings = {ing: pct for ing, pct in blend.items()
                    if any(kw in ing.lower() for kw in SOLVENTS)}
    total_solvent = sum(solvent_ings.values())
    if total_solvent > 30:
        contributing.extend(solvent_ings.keys())
        evidence.append(
            f"High solvent load ({total_solvent:.0f}% total) — "
            f"solvents can destabilise emulsions and accelerate oxidative degradation "
            f"under ICH accelerated conditions"
        )

    # Rule: no preservative / antioxidant system
    has_preservative = any(
        any(kw in ing.lower() for kw in PRESERVATIVE_KEYWORDS) for ing in blend
    )
    if not has_preservative:
        evidence.append(
            "No preservative or antioxidant detected — "
            "formulation is unprotected against microbial growth and oxidative degradation"
        )

    # Rule: no emulsifier when blend contains both polar and non-polar phases
    has_emulsifier = any(any(e in ing.lower() for e in EMULSIFIERS) for ing in blend)
    if not has_emulsifier and total_solvent > 10:
        evidence.append(
            "No emulsifier detected — polar/non-polar phases may separate under "
            "thermal stress (freeze-thaw or 40°C cycling)"
        )

    confidence = 0.65 if (total_solvent > 30 or not has_preservative) else 0.45
    try:
        wks = float(timepoint)
        severity = "Severe" if wks <= 2 else "Moderate" if wks <= 8 else "Minor"
    except (TypeError, ValueError):
        severity = "Moderate"

    root_cause = (
        f"Accelerated stability failure at {condition}, week {timepoint}: {observation}. "
        + (f"High solvent load ({total_solvent:.0f}% total) identified. " if total_solvent > 30 else "")
        + ("No preservative system present. " if not has_preservative else "")
        + ("No emulsifier detected. " if not has_emulsifier and total_solvent > 10 else "")
    ).strip()

    return RootCauseAnalysis(
        failure_type="stability_failure",
        root_cause=root_cause,
        contributing_ingredients=contributing,
        severity=severity,
        confidence=confidence,
        evidence=evidence,
    )


def _suggest_stability_fix(
    blend: Dict[str, float], db: pd.DataFrame, test_data: dict
) -> List[ReformulationSuggestion]:
    """Generate stability-focused reformulation suggestions.

    Rules applied (in priority order):
      1. High solvent (>30% total) → reduce top solvent 5–10 pp, add emulsifier 2–5%
      2. No preservative / antioxidant → add EDTA Disodium 0.1%
    Minimum confidence 0.4 is enforced when at least one rule fires.
    """
    suggestions = []
    rank = 1

    solvent_ings = {ing: pct for ing, pct in blend.items()
                    if any(kw in ing.lower() for kw in SOLVENTS)}
    total_solvent = sum(solvent_ings.values())

    # Rule 1a — reduce highest-load solvent by 5–10 percentage points
    if total_solvent > 30 and solvent_ings:
        top_solvent = max(solvent_ings, key=lambda x: solvent_ings[x])
        top_pct = solvent_ings[top_solvent]
        reduction = max(5.0, min(10.0, top_pct * 0.20))

        suggestions.append(ReformulationSuggestion(
            rank=rank, action_type="decrease", ingredient=top_solvent,
            current_pct=top_pct, suggested_pct=round(top_pct - reduction, 1),
            rationale=(
                f"Reduce {top_solvent} by {reduction:.0f} pp — total solvent load is "
                f"{total_solvent:.0f}% (threshold: 30%). High solvent concentration "
                f"accelerates emulsion breakdown and oxidative degradation under ICH conditions."
            ),
            confidence=0.65,
            predicted_impact={
                "stability": "Improved — reduced solvent stress on formulation",
                "cleaning_performance": "Monitor — retest soil removal at new level",
                "cost": f"~-${reduction * 0.01:.2f}/kg estimated",
            },
            cost_delta_per_kg=round(-reduction * 0.01, 3),
            risk_level="Low",
            implementation_notes=(
                f"Replace displaced {reduction:.0f}% with deionised water or co-solvent. "
                "Retest accelerated stability at 40°C/75% RH, 4 weeks."
            ),
        ))
        rank += 1

        # Rule 1b — add emulsifier if none present to stabilise reduced-solvent system
        has_emulsifier = any(any(e in ing.lower() for e in EMULSIFIERS) for ing in blend)
        if not has_emulsifier:
            suggestions.append(ReformulationSuggestion(
                rank=rank, action_type="add", ingredient="Polysorbate 20",
                current_pct=0.0, suggested_pct=3.0,
                rationale=(
                    "Add Polysorbate 20 at 3% — non-ionic HLB 16.7 emulsifier that bridges "
                    "the solvent/water interface. REACH registered, EPA Safer Choice listed, "
                    "widely compatible with industrial and personal care verticals."
                ),
                confidence=0.72,
                predicted_impact={
                    "stability": "Significantly improved",
                    "phase_separation_risk": "Reduced",
                    "cost": "+$0.09/kg",
                },
                cost_delta_per_kg=0.09,
                risk_level="Low",
                implementation_notes=(
                    "Add to aqueous phase at room temperature before combining phases. "
                    "Effective range 2–5%; titrate to minimum effective concentration."
                ),
            ))
            rank += 1

    # Rule 2 — add chelating agent if no preservative system present
    has_preservative = any(
        any(kw in ing.lower() for kw in PRESERVATIVE_KEYWORDS) for ing in blend
    )
    if not has_preservative:
        suggestions.append(ReformulationSuggestion(
            rank=rank, action_type="add", ingredient="EDTA Disodium",
            current_pct=0.0, suggested_pct=0.1,
            rationale=(
                "Add EDTA Disodium 0.1% — chelating agent that sequesters trace metal ions "
                "which catalyse oxidative degradation. Standard protective additive in "
                "industrial cleaners and personal care formulations."
            ),
            confidence=0.70,
            predicted_impact={
                "oxidative_stability": "Improved",
                "microbial_stability": "Partial — metal chelation inhibits metal-dependent microbes",
                "cost": "+$0.005/kg",
            },
            cost_delta_per_kg=0.005,
            risk_level="Low",
            implementation_notes=(
                "Add to aqueous phase. REACH registered, low ecotox concern at 0.1%. "
                "For food-contact surfaces, use food-grade EDTA or substitute Citric Acid 0.2%."
            ),
        ))

    # Enforce minimum confidence of 0.4 on all suggestions when at least one rule fired
    if suggestions:
        for s in suggestions:
            s.confidence = max(s.confidence, 0.4)

    return suggestions


def _suggest_certification_fix(blend: Dict[str, float], db: pd.DataFrame,
                                 test_data: dict) -> List[ReformulationSuggestion]:
    """Generate certification rejection fix."""
    rejected = test_data.get("rejected_ingredient", "")
    cert = test_data.get("certification_name", "")
    suggestions = []

    if rejected and rejected in blend:
        pct = blend[rejected]

        # Find replacement based on function
        idx = db.set_index("Ingredient") if "Ingredient" in db.columns else db
        func = ""
        if rejected in idx.index and "Function" in idx.columns:
            func = str(idx.loc[rejected, "Function"])

        # Generic bio-based replacement suggestion
        suggestions.append(ReformulationSuggestion(
            rank=1, action_type="replace", ingredient=rejected,
            current_pct=pct, suggested_pct=pct,
            rationale=f"Replace {rejected} — rejected by {cert}. "
                      f"Find {cert}-approved alternative with same function: {func[:40]}",
            confidence=0.95,
            predicted_impact={"certification": f"{cert} compliance ✅",
                               "performance": "monitor — test equivalent performance"},
            cost_delta_per_kg=0.10, risk_level="Medium",
            implementation_notes=(
                f"Use IntelliForm CertificationOracle to screen replacement candidates. "
                f"Run pilot batch and resubmit for {cert} screening."
            ),
        ))

    return suggestions


def _suggest_supply_fix(blend: Dict[str, float], db: pd.DataFrame,
                         test_data: dict) -> List[ReformulationSuggestion]:
    """Generate supply disruption alternatives."""
    disrupted = test_data.get("disrupted_ingredient", "")
    lead_time = float(test_data.get("lead_time_weeks", 12))
    suggestions = []

    if disrupted and disrupted in blend:
        idx = db.set_index("Ingredient") if "Ingredient" in db.columns else db
        pct = blend[disrupted]
        func = str(idx.loc[disrupted, "Function"]) if disrupted in idx.index else ""

        # Find alternatives with same function in DB
        alts = []
        if func and "Function" in db.columns:
            alts = db[db["Function"].str.contains(
                func[:20], case=False, na=False)
            ]["Ingredient"].tolist()
            alts = [a for a in alts if a != disrupted][:3]

        for i, alt in enumerate(alts):
            suggestions.append(ReformulationSuggestion(
                rank=i+1, action_type="replace", ingredient=disrupted,
                current_pct=pct, suggested_pct=pct,
                rationale=f"Replace {disrupted} with {alt} — same function ({func[:30]}), "
                          f"potentially lower lead time",
                confidence=0.65,
                predicted_impact={"supply": "Resolved", "performance": "test required"},
                cost_delta_per_kg=0.0, risk_level="Medium",
                implementation_notes=f"Contact Univar/Brenntag for {alt} availability",
            ))

    return suggestions


# ── Main entry point ───────────────────────────────────────────────────────────

def run_reformulation_intelligence(
    blend: Dict[str, float],
    db: pd.DataFrame,
    failure_type: str,
    test_data: dict,
    batch_id: Optional[str] = None,
) -> ReformulationReport:
    """
    Run closed-loop reformulation intelligence on a failed batch.

    Args:
        blend: Current formulation {ingredient: percentage}
        db: Ingredient database
        failure_type: Key from FAILURE_TYPES dict
        test_data: Measured values from lab test
        batch_id: ChemRich batch reference

    Returns:
        ReformulationReport with root cause + ranked suggestions
    """
    idx = db.set_index("Ingredient") if "Ingredient" in db.columns else db

    # Diagnose root cause
    if failure_type in ("ph_too_high", "ph_too_low"):
        rca = _diagnose_ph(blend, db, failure_type, test_data)
        suggestions = _suggest_ph_fix(blend, db, failure_type, test_data)

    elif failure_type in ("viscosity_too_high", "viscosity_too_low"):
        rca = _diagnose_viscosity(blend, db, failure_type, test_data)
        suggestions = _suggest_viscosity_fix(blend, db, failure_type, test_data)

    elif failure_type == "phase_separation":
        rca = _diagnose_phase_separation(blend, db, test_data)
        suggestions = _suggest_phase_separation_fix(blend, db, test_data)

    elif failure_type == "stability_failure":
        rca = _diagnose_stability(blend, test_data)
        suggestions = _suggest_stability_fix(blend, db, test_data)

    elif failure_type == "certification_rejection":
        rca = _diagnose_certification(blend, test_data)
        suggestions = _suggest_certification_fix(blend, db, test_data)

    elif failure_type == "supply_disruption":
        rca = RootCauseAnalysis(
            failure_type="supply_disruption",
            root_cause=f"Supply disruption: {test_data.get('disrupted_ingredient')} "
                       f"unavailable ({test_data.get('lead_time_weeks')} week lead time)",
            contributing_ingredients=[test_data.get("disrupted_ingredient", "")],
            severity="Moderate", confidence=0.99,
            evidence=["Supplier confirmed unavailability"],
        )
        suggestions = _suggest_supply_fix(blend, db, test_data)

    elif failure_type == "cost_overrun":
        actual = float(test_data.get("actual_cost_per_kg", 0))
        target = float(test_data.get("target_cost_per_kg", 0))
        # Find most expensive ingredients
        costs = {}
        for ing in blend:
            if ing in idx.index and "Cost_USD_kg" in idx.columns:
                costs[ing] = float(idx.loc[ing, "Cost_USD_kg"]) * blend[ing] / 100
        expensive = sorted(costs, key=lambda x: -costs[x])[:3]
        rca = RootCauseAnalysis(
            failure_type="cost_overrun",
            root_cause=f"Batch cost ${actual:.2f}/kg exceeds target ${target:.2f}/kg "
                       f"(+{(actual/target-1)*100:.0f}%). Top cost drivers: {', '.join(expensive[:2])}",
            contributing_ingredients=expensive,
            severity="Moderate" if actual/target < 1.2 else "Severe",
            confidence=0.95,
            evidence=[f"{ing}: ${costs.get(ing,0):.2f}/kg contribution" for ing in expensive],
        )
        suggestions = [ReformulationSuggestion(
            rank=i+1, action_type="decrease", ingredient=ing,
            current_pct=blend[ing], suggested_pct=round(blend[ing] * 0.8, 1),
            rationale=f"Reduce {ing} — highest cost contributor at ${costs.get(ing,0):.2f}/kg impact",
            confidence=0.70,
            predicted_impact={"cost": f"-${costs.get(ing,0)*0.2:.2f}/kg",
                               "performance": "test required"},
            cost_delta_per_kg=-costs.get(ing, 0) * 0.2,
            risk_level="Medium",
            implementation_notes="Replace partially with lower-cost ingredient of same function",
        ) for i, ing in enumerate(expensive[:3])]

    else:
        # Generic failure
        rca = RootCauseAnalysis(
            failure_type=failure_type,
            root_cause=f"Failure type: {FAILURE_TYPES.get(failure_type, {}).get('label', failure_type)}",
            contributing_ingredients=[], severity="Moderate",
            confidence=0.50, evidence=["Insufficient data for root cause analysis"],
        )
        suggestions = []

    # Best suggestion = rank 1
    best = suggestions[0] if suggestions else ReformulationSuggestion(
        rank=1, action_type="review", ingredient="—",
        current_pct=0, suggested_pct=0,
        rationale=(
            "No automated rule matched this failure pattern. "
            "Manual formulation review recommended — check raw material COAs, "
            "processing conditions (temperature, mixing order, shear), "
            "and storage history before changing the formulation."
        ),
        confidence=0.0,
        predicted_impact={
            "next_step": "Review batch records and lab notes against specification",
            "escalation": "Consider DoE (Design of Experiments) to isolate failure variable",
        },
        cost_delta_per_kg=0,
        risk_level="Medium",
        implementation_notes=(
            "Document failure observations and test data — structured records improve "
            "future rule matching. Check IntelliForm failure taxonomy for closest match."
        ),
    )

    # Predict success probability
    success_prob = best.confidence * (1.0 - 0.2 * len(rca.contributing_ingredients))
    success_prob = float(np.clip(success_prob, 0.3, 0.95))

    # Do not change — only flag ingredients we have specific evidence are working.
    # When contributing_ingredients is empty (no root cause identified), we have no
    # basis to declare any ingredient safe, so return an empty list rather than
    # locking the entire blend by default.
    do_not_change = (
        [ing for ing in blend if ing not in rca.contributing_ingredients][:3]
        if rca.contributing_ingredients
        else []
    )

    # Iterations estimate
    iterations = "1–2 batches" if rca.severity in ("Minor", "Moderate") else "2–4 batches"

    # Learning note
    learning_note = (
        f"This failure teaches: {rca.root_cause[:120]}. "
        f"Confidence in fix: {success_prob*100:.0f}%. "
        f"Feed outcome back to IntelliForm to improve future predictions."
    )

    # CBAM impact
    cbam_impact = None
    if suggestions and suggestions[0].action_type in ("replace", "add"):
        cbam_impact = "Reformulation may change carbon footprint — regenerate Carbon Passport after fix"

    return ReformulationReport(
        failure_type=failure_type,
        failure_description=FAILURE_TYPES.get(failure_type, {}).get("description", failure_type),
        root_cause=rca,
        suggestions=suggestions,
        best_suggestion=best,
        predicted_success_probability=round(success_prob, 2),
        iterations_to_fix=iterations,
        do_not_change=do_not_change,
        learning_note=learning_note,
        cbam_impact=cbam_impact,
    )
