"""
modules/stability.py
Stability and viscosity prediction for IntelliForm v1.0.

Predicts:
  - Shelf life (months) based on blend composition and preservation system
  - Viscosity (cP at 25C) from thickener content and surfactant ratios
  - pH stability range
  - Recommended packaging type

Rule-based models calibrated against published formulation literature.
Upgrades to ML when lab validation data is available.
"""
from dataclasses import dataclass, field
from typing import Dict, List, Optional
import pandas as pd


@dataclass
class StabilityResult:
    shelf_life_months: int           # predicted shelf life
    shelf_life_range: str            # e.g. "18-24 months"
    viscosity_cp: float              # predicted viscosity at 25C
    viscosity_range: str             # e.g. "2,800-3,600 cP"
    ph_min: float                    # stable pH lower bound
    ph_max: float                    # stable pH upper bound
    recommended_packaging: str       # HDPE, glass, PET, aluminium
    stability_risks: List[str]       # flagged risks
    stability_boosters: List[str]    # positive stability factors
    overall_rating: str              # Excellent / Good / Fair / Poor
    confidence: str                  # high / medium / low


# ── Function-based lookup tables ──────────────────────────────────────────────

_THICKENER_VISCOSITY = {
    "Xanthan Gum": 8000,
    "Hydroxyethyl Cellulose": 4000,
    "Guar Gum (bio)": 5000,
    "Hydroxypropyl Methylcellulose": 3500,
    "Carrageenan": 3000,
    "Cellulose Gum (CMC)": 2500,
    "Pectin": 2000,
    "Carbomer": 9000,
    "Sodium Alginate": 3000,
    "Hydroxypropyl Starch": 2000,
}

_PRESERVATIVE_EFFICACY = {
    "Benzyl Alcohol (NF)": 18,
    "Sodium Benzoate": 12,
    "Potassium Sorbate": 12,
    "Phenoxyethanol": 24,
    "Ethylhexylglycerin": 18,
    "Caprylhydroxamic Acid": 24,
    "Caprylic Acid (bio)": 18,
    "Undecylenic Acid": 18,
    "Levulinic Acid (bio)": 15,
    "Hydroxyacetophenone": 18,
    "1-Hexanediol": 18,
    "Hexylene Glycol": 18,
}

_ANTIOXIDANT_BOOST = {
    "Tocopheryl Acetate": 6,
    "Ascorbic Acid": 4,
    "Vitamin E (mixed tocopherols)": 8,
    "Rosemary Extract (antioxidant)": 6,
    "Green Tea Extract": 5,
}

_EMOLLIENT_PACKAGING = {
    "Argan Oil": "Glass or aluminium — protect from UV",
    "Jojoba Esters": "HDPE or glass acceptable",
    "Coconut Oil (fractionated)": "HDPE",
    "Shea Butter (refined)": "HDPE or PET",
}

_RISK_INGREDIENTS = {
    "D-Limonene": "D-Limonene can oxidise on prolonged storage — add antioxidant",
    "Rhamnolipid (biosurfactant)": "Biosurfactants sensitive to pH shift — monitor monthly",
    "Sophorolipid": "Sophorolipids may crystallise below 15C — store above ambient",
    "Sodium Coco-Sulfate": "High sulfate content may cause pH drift over time",
    "Lecithin (sunflower)": "Lecithin susceptible to hydrolysis — keep pH 4-7",
    "Lecithin (soy)": "Lecithin susceptible to hydrolysis — keep pH 4-7",
    "Ascorbic Acid": "Ascorbic Acid oxidises rapidly — use nitrogen headspace packaging",
    "Ceramide NP (bio)": "Ceramides require emulsification stabiliser for aqueous systems",
}


# ── Core prediction ───────────────────────────────────────────────────────────

def predict_stability(blend: Dict[str, float], db: pd.DataFrame) -> StabilityResult:
    """
    Predict stability and viscosity for a formulation blend.
    """
    idx = db.set_index("Ingredient") if "Ingredient" in db.columns else db

    # ── Viscosity prediction ──
    base_viscosity = 50.0  # water-like base
    for ing, pct in blend.items():
        if ing in _THICKENER_VISCOSITY:
            # Viscosity contribution scales non-linearly with concentration
            conc_factor = (pct / 100) ** 0.7
            base_viscosity += _THICKENER_VISCOSITY[ing] * conc_factor

    # Surfactants reduce viscosity slightly
    for ing, pct in blend.items():
        try:
            if ing in idx.index:
                func = str(idx.loc[ing, "Function"]) if "Function" in idx.columns else ""
                if "Surfactant" in func:
                    base_viscosity *= (1 - (pct / 100) * 0.15)
        except Exception:
            pass

    viscosity = round(base_viscosity, 0)
    v_low = round(viscosity * 0.85, 0)
    v_high = round(viscosity * 1.15, 0)

    # ── Shelf life prediction ──
    base_shelf_life = 12  # months — conservative default

    # Preservative system
    max_preservative = 0
    for ing, pct in blend.items():
        if ing in _PRESERVATIVE_EFFICACY and pct > 0.5:
            max_preservative = max(max_preservative, _PRESERVATIVE_EFFICACY[ing])
    base_shelf_life = max(base_shelf_life, max_preservative)

    # Antioxidant boost
    antioxidant_boost = 0
    for ing in blend:
        if ing in _ANTIOXIDANT_BOOST:
            antioxidant_boost = max(antioxidant_boost, _ANTIOXIDANT_BOOST[ing])
    base_shelf_life += antioxidant_boost

    # High bio-based content slightly reduces shelf life (more biodegradable)
    try:
        bio_avg = sum(
            float(idx.loc[ing, "Bio_based_pct"]) * (pct / 100)
            for ing, pct in blend.items() if ing in idx.index
        )
        if bio_avg > 95:
            base_shelf_life = max(base_shelf_life - 3, 12)
    except Exception:
        pass

    shelf_low = max(base_shelf_life - 3, 6)
    shelf_high = base_shelf_life + 3

    # ── pH range ──
    ph_min, ph_max = 4.5, 8.5
    for ing in blend:
        if "Citric Acid" in ing or "Lactic Acid" in ing or "Glucono" in ing:
            ph_min, ph_max = 4.0, 6.5
        if "Sodium Metasilicate" in ing or "Sodium Carbonate" in ing:
            ph_min, ph_max = 9.0, 12.0
        if "Sodium Bicarbonate" in ing:
            ph_min, ph_max = 7.5, 9.0

    # ── Packaging recommendation ──
    packaging = "HDPE (recommended for most green chemistry formulations)"
    for ing in blend:
        if ing in _EMOLLIENT_PACKAGING:
            packaging = _EMOLLIENT_PACKAGING[ing]
            break
    if ph_min >= 9.0:
        packaging = "HDPE only — alkaline formulas degrade PET and glass closures"
    if any("Ascorbic Acid" in ing or "Vitamin C" in ing for ing in blend):
        packaging = "Aluminium or dark glass — protect from UV and oxygen"

    # ── Risk flags ──
    risks = []
    for ing in blend:
        if ing in _RISK_INGREDIENTS:
            risks.append(_RISK_INGREDIENTS[ing])
    if not risks:
        risks = ["No significant stability risks identified"]

    # ── Stability boosters ──
    boosters = []
    if antioxidant_boost > 0:
        boosters.append("Antioxidant system detected — extends oxidative stability")
    if max_preservative >= 18:
        boosters.append("Broad-spectrum preservative system — good microbial protection")
    for ing in blend:
        if "Chelating" in str(idx.loc[ing, "Function"]) if ing in idx.index and "Function" in idx.columns else "":
            boosters.append(f"{ing} chelates heavy metals — prevents catalytic oxidation")
            break
    if not boosters:
        boosters = ["Consider adding antioxidant and chelating agent to extend shelf life"]

    # ── Overall rating ──
    score = 0
    if base_shelf_life >= 24: score += 3
    elif base_shelf_life >= 18: score += 2
    else: score += 1
    if max_preservative > 0: score += 2
    if antioxidant_boost > 0: score += 1
    if len(risks) <= 1: score += 1

    if score >= 6: rating = "Excellent"
    elif score >= 4: rating = "Good"
    elif score >= 2: rating = "Fair"
    else: rating = "Poor"

    return StabilityResult(
        shelf_life_months=base_shelf_life,
        shelf_life_range=f"{shelf_low}–{shelf_high} months",
        viscosity_cp=viscosity,
        viscosity_range=f"{v_low:,.0f}–{v_high:,.0f} cP",
        ph_min=ph_min,
        ph_max=ph_max,
        recommended_packaging=packaging,
        stability_risks=risks,
        stability_boosters=boosters,
        overall_rating=rating,
        confidence="medium",
    )
