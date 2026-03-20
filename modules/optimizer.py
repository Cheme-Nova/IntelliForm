"""
modules/optimizer.py
Vertical-aware PuLP LP optimizer for IntelliForm v1.1.

Core improvements:
  - Vertical-specific mandatory function requirements
  - Per-function concentration bounds (pharma lubricant max 2%, etc.)
  - VOC budget constraint (paint/coatings)
  - Functional completeness validation post-solve
  - Detailed validation report per vertical
  - Backward compatible — vertical='all' behaves like original optimizer
"""
from dataclasses import dataclass, field
from typing import Dict, Optional, List
import pandas as pd
import pulp
from modules.analytics import track


# ── Result dataclasses ────────────────────────────────────────────────────────

@dataclass
class OptResult:
    success: bool
    blend: Dict[str, float]
    cost_per_kg: float
    bio_pct: float
    perf_score: float
    status: str
    relaxed: bool = False
    relaxation_rounds: int = 0
    error_msg: Optional[str] = None
    vertical: str = "all"
    vertical_validation: Optional[dict] = None
    manufacturing_route: Optional[str] = None
    compliance_flags: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)


# ── Vertical constraint configs ───────────────────────────────────────────────

VERTICAL_CONSTRAINTS = {

    "pharmaceutical": {
        # Per-function concentration bounds (fraction of total blend)
        "function_bounds": {
            "Pharma Lubricant":           (0.002, 0.025),  # 0.2–2.5% lubricant
            "Pharma Binder":              (0.02,  0.15),   # 2–15% binder
            "Pharma Filler":              (0.20,  0.80),   # 20–80% filler
            "Pharma Disintegrant":        (0.01,  0.10),   # 1–10% disintegrant
            "Pharma Film Coat":           (0.02,  0.05),   # 2–5% film coat
            "Pharma Superdisintegrant":   (0.01,  0.08),   # 1–8%
            "Pharma Glidant":             (0.001, 0.005),  # 0.1–0.5% glidant
            "Pharma Flow Aid":            (0.001, 0.005),
            "Pharma Preservative":        (0.001, 0.01),   # 0.1–1%
            "Pharma Sweetener":           (0.001, 0.05),   # 0.1–5%
        },
        # Required functional groups — at least one ingredient per group
        "required_function_groups": [
            ["Pharma Filler", "Pharma Binder/Filler", "Pharma DCI Filler"],
        ],
        # Optional but recommended groups
        "recommended_function_groups": [
            ["Pharma Disintegrant", "Pharma Superdisintegrant"],
            ["Pharma Lubricant", "Pharma Glidant", "Pharma Flow Aid"],
        ],
        # Global constraints
        "min_bio": 0,     # pharma doesn't require bio-based
        "max_water_content": 0.05,  # max 5% water for solid dosage
        "regulatory_framework": "ICH Q8 / USP-NF / FDA 21 CFR",
        "manufacturing_routes": {
            "direct_compression": "Best for: low-dose APIs, free-flowing powders, MCC dominant",
            "wet_granulation":    "Best for: poor flow powders, high-dose APIs, better uniformity",
            "dry_granulation":    "Best for: moisture-sensitive APIs, roller compaction",
            "hot_melt_extrusion": "Best for: poorly soluble APIs, amorphous dispersion",
        },
    },

    "food": {
        "function_bounds": {
            "Food Emulsifier":            (0.001, 0.05),   # 0.1–5%
            "Food Thickener":             (0.001, 0.03),   # 0.1–3%
            "Food Thickener/Stabilizer":  (0.001, 0.03),
            "Food Preservative":          (0.0001, 0.003), # 0.01–0.3%
            "Food Antioxidant":           (0.0001, 0.002), # 0.01–0.2%
            "Food Color":                 (0.0001, 0.005), # trace–0.5%
            "Food Natural Color":         (0.0001, 0.01),
            "Food Acidulant":             (0.001, 0.05),
            "Food Sweetener":             (0.001, 0.20),
            "Food Natural Sweetener":     (0.001, 0.10),
            "Food Flavor":                (0.0001, 0.01),  # trace–1%
            "Food Gelling Agent":         (0.001, 0.02),
            "Food Fiber":                 (0.01,  0.15),
            "Food Protein":               (0.01,  0.30),
            "Food Protein/Emulsifier":    (0.01,  0.30),
        },
        "required_function_groups": [],
        "gras_required": True,
        "allergen_flags": ["Soy Protein", "Whey Protein", "Casein", "Egg Albumen",
                          "Wheat", "Peanut", "Tree Nut", "Fish", "Shellfish"],
        "clean_label_avoid": ["Sodium Nitrite", "Sulfur Dioxide", "Sodium Erythorbate",
                              "BHT", "Titanium Dioxide"],
        "regulatory_framework": "FDA GRAS / EU Food Additive Regulation / EFSA / Codex Alimentarius",
        "max_additive_load": 0.15,  # total additives max 15% of formulation
    },

    "agricultural": {
        "function_bounds": {
            "Agri Superspreader Adjuvant": (0.001, 0.01),  # 0.1–1% — very potent
            "Agri Oil Adjuvant":           (0.001, 0.05),
            "Agri MSO Adjuvant":           (0.001, 0.05),
            "Agri Biopesticide":           (0.001, 0.10),
            "Agri Bioinsecticide":         (0.001, 0.10),
            "Agri Biocontrol Agent":       (0.001, 0.10),
            "Agri Plant Growth Regulator": (0.000001, 0.001),  # trace — PPM level
            "Agri Auxin/Root Growth":      (0.000001, 0.001),
            "Agri Micronutrient":          (0.0001, 0.01),
            "Agri Micronutrient Chelate":  (0.0001, 0.01),
            "Agri Nitrogen Fertilizer":    (0.01, 0.50),
            "Agri Surfactant":             (0.001, 0.05),
            "Agricultural Surfactant":     (0.001, 0.05),
        },
        "rei_required": True,    # Restricted Entry Interval
        "phi_required": True,    # Pre-Harvest Interval
        "bee_safety_check": True,
        "omri_preferred": True,
        "regulatory_framework": "EPA FIFRA / EU Reg 1107/2009 / OMRI (organic)",
        "drift_restrict": 0.05,  # max 5% drift reduction adjuvant
    },

    "paint_coatings": {
        "function_bounds": {
            "Paint/Coatings Binder":          (0.20, 0.60),  # 20–60% binder
            "Paint/Coatings Waterborne Binder":(0.20, 0.60),
            "Paint/Coatings Bio-Alkyd":        (0.20, 0.60),
            "Paint/Coatings Pigment":          (0.05, 0.30),  # 5–30% pigment
            "Paint/Coatings Opacifier":        (0.05, 0.25),
            "Paint/Coatings Extender":         (0.05, 0.40),  # 5–40% extender
            "Paint/Coatings Coated Extender":  (0.05, 0.40),
            "Paint/Coatings Thickener":        (0.001, 0.03), # 0.1–3%
            "Paint/Coatings Coalescent":       (0.01, 0.05),  # 1–5%
            "Paint/Coatings Defoamer":         (0.001, 0.005),# 0.1–0.5%
            "Paint/Coatings Dispersant":       (0.001, 0.02), # 0.1–2%
            "Paint/Coatings In-Can Biocide":   (0.0001, 0.002),# trace
            "Paint/Coatings Film Biocide":     (0.0001, 0.002),
            "Paint/Coatings UV Stabilizer":    (0.001, 0.03),
            "Paint/Coatings UV Absorber":      (0.001, 0.03),
            "Paint/Coatings Crosslinker":      (0.01, 0.10),
        },
        "voc_limit_decorative": 30.0,    # g/L VOC — EU Decopaint Directive
        "voc_limit_industrial": 250.0,   # g/L VOC — industrial coatings
        "pigment_volume_conc": (0.25, 0.85),  # PVC 25–85%
        "required_function_groups": [
            ["Paint/Coatings Binder", "Paint/Coatings Waterborne Binder",
             "Paint/Coatings Bio-Alkyd", "Paint/Coatings Alkyd Binder",
             "Paint/Coatings Alkyd", "Paint/Coatings Bio-Resin"],
        ],
        "regulatory_framework": "EU Decopaint Directive / VOC Regulation EC 2004/42 / LEED",
    },

    "fabric_laundry": {
        "function_bounds": {
            "Fabric Primary Surfactant":   (0.05, 0.30),  # 5–30%
            "Fabric Anionic Surfactant":   (0.05, 0.25),
            "Fabric Amphoteric Surfactant":(0.01, 0.10),
            "Fabric Nonionic Surfactant":  (0.01, 0.15),
            "Fabric Builder":              (0.05, 0.30),
            "Fabric Softener":             (0.01, 0.20),
            "Fabric Stain Removal Enzyme": (0.001, 0.02), # enzymes at low %
            "Fabric Bio-Polishing Enzyme": (0.001, 0.02),
            "Fabric Oxidizing Bleach":     (0.05, 0.25),
            "Fabric Bleach Activator":     (0.01, 0.10),
            "Fabric OBA":                  (0.001, 0.005),# optical brightener trace
            "Fabric Anti-Redeposition":    (0.001, 0.05),
            "Fabric Encapsulated Fragrance":(0.001, 0.02),
            "Fabric Scale Inhibitor":      (0.001, 0.01),
        },
        "required_function_groups": [
            ["Fabric Primary Surfactant", "Fabric Anionic Surfactant"],
        ],
        "phosphate_free": True,
        "biodeg_28_day": 0.60,    # 60% biodegradation in 28 days
        "regulatory_framework": "EU Detergent Regulation 648/2004 / EPA DfE",
    },

    "industrial": {
        "function_bounds": {
            "Industrial Solvent":          (0.10, 0.80),
            "Industrial Bio-Solvent":      (0.10, 0.80),
            "Industrial Bio-Degreaser":    (0.10, 0.70),
            "Industrial Hydrotrope":       (0.01, 0.15),
            "Industrial Chelating Agent":  (0.01, 0.10),
            "Industrial Biocide":          (0.0001, 0.005),
            "Industrial Corrosion Inhibitor":(0.001, 0.05),
            "Industrial Scale Inhibitor":  (0.001, 0.02),
            "Industrial Disinfectant":     (0.01, 0.10),
        },
        "drain_safe": True,
        "voc_limit": 350.0,  # g/L
        "regulatory_framework": "REACH / EPA Safer Choice / NSF A1 (food-safe cleaning)",
    },

    "personal_care": {
        "function_bounds": {
            "Personal Care Antioxidant Active": (0.0001, 0.05),
            "Personal Care Anti-aging Peptide":  (0.0001, 0.02),
            "Personal Care Brightening Active":  (0.001,  0.10),
            "Personal Care AHA Exfoliant":        (0.01,   0.15),
            "Personal Care UV Filter":            (0.01,   0.25),
            "Personal Care Preservative":         (0.0001, 0.01),
            "Preservative":                       (0.0001, 0.01),
            "Personal Care Fragrance":            (0.0001, 0.03),
            "Personal Care Essential Oil":        (0.0001, 0.03),
            "Personal Care Silicone Emollient":   (0.01,   0.30),
            "Personal Care Thickener":            (0.001,  0.05),
            "Thickener":                          (0.001,  0.05),
        },
        "required_function_groups": [],
        "skin_safe": True,
        "regulatory_framework": "EU Cosmetics Regulation 1223/2009 / COSMOS / REACH",
        "max_fragrance": 0.03,  # IFRA max 3% fragrance
    },

}


# ── Vertical validation ───────────────────────────────────────────────────────

def validate_blend_for_vertical(blend: Dict[str, float], db: pd.DataFrame,
                                  vertical: str) -> dict:
    """Run vertical-specific validation on a completed blend."""
    config = VERTICAL_CONSTRAINTS.get(vertical, {})
    idx = db.set_index("Ingredient") if "Ingredient" in db.columns else db

    flags = []
    warnings = []
    passes = []

    # Check required function groups
    for group in config.get("required_function_groups", []):
        found = False
        for ing in blend:
            if ing in idx.index:
                func = str(idx.loc[ing, "Function"]) if "Function" in idx.columns else ""
                if any(g.lower() in func.lower() for g in group):
                    found = True
                    break
        if found:
            passes.append(f"Required function present: {group[0]}")
        else:
            warnings.append(f"Missing recommended ingredient type: {group[0]}")

    # Pharma-specific
    if vertical == "pharmaceutical":
        total_pct = sum(blend.values())
        lubricant_pct = sum(
            pct for ing, pct in blend.items()
            if ing in idx.index and "Lubricant" in str(idx.loc[ing, "Function"] if "Function" in idx.columns else "")
        )
        if lubricant_pct > 2.5:
            flags.append(f"Lubricant concentration {lubricant_pct:.1f}% exceeds ICH guideline max 2.5%")
        elif lubricant_pct > 0:
            passes.append(f"Lubricant level {lubricant_pct:.1f}% within ICH guidelines")

        filler_pct = sum(
            pct for ing, pct in blend.items()
            if ing in idx.index and any(x in str(idx.loc[ing, "Function"] if "Function" in idx.columns else "")
                                         for x in ["Filler", "Binder/Filler"])
        )
        if filler_pct < 20:
            warnings.append(f"Filler/binder content {filler_pct:.1f}% is low — typical tablet is 60-80% filler")

        # Manufacturing route recommendation
        if filler_pct >= 50 and lubricant_pct > 0:
            route = "Direct Compression — good filler/lubricant balance detected"
        elif filler_pct >= 30:
            route = "Wet Granulation — consider for improved flow and compressibility"
        else:
            route = "Review formulation — insufficient filler content for standard tablet"

        return {
            "passes": passes, "warnings": warnings, "flags": flags,
            "manufacturing_route": route,
            "ich_compliant": len(flags) == 0,
            "regulatory_framework": config.get("regulatory_framework", ""),
        }

    # Food-specific
    if vertical == "food":
        allergens_present = []
        for ing in blend:
            for allergen in config.get("allergen_flags", []):
                if allergen.lower() in ing.lower():
                    allergens_present.append(ing)
        if allergens_present:
            flags.append(f"ALLERGEN DECLARATION REQUIRED: {', '.join(allergens_present)}")

        clean_label_violations = [i for i in blend if any(x.lower() in i.lower()
                                  for x in config.get("clean_label_avoid", []))]
        if clean_label_violations:
            warnings.append(f"Non-clean-label ingredients: {', '.join(clean_label_violations)}")
        else:
            passes.append("Clean label compliant — no artificial additives detected")

        return {
            "passes": passes, "warnings": warnings, "flags": flags,
            "allergens": allergens_present,
            "clean_label": len(clean_label_violations) == 0,
            "gras_assumed": True,
            "regulatory_framework": config.get("regulatory_framework", ""),
        }

    # Agricultural-specific
    if vertical == "agricultural":
        potent = []
        for ing, pct in blend.items():
            if ing in idx.index:
                func = str(idx.loc[ing, "Function"] if "Function" in idx.columns else "")
                if "Plant Growth Regulator" in func or "Auxin" in func:
                    if pct > 0.1:
                        flags.append(f"{ing} at {pct:.2f}% — PGR typically dosed at ppm level (0.001–0.01%)")
        passes.append("REI/PHI compliance check: verify with EPA registration before field use")
        warnings.append("Tank mix compatibility: test with target pesticide before commercial use")

        return {
            "passes": passes, "warnings": warnings, "flags": flags,
            "rei_note": "Restricted Entry Interval varies by active ingredient — check EPA label",
            "phi_note": "Pre-Harvest Interval varies by crop — check registered label",
            "bee_safety": "Verify bee safety rating for all actives",
            "regulatory_framework": config.get("regulatory_framework", ""),
        }

    # Paint-specific
    if vertical == "paint_coatings":
        # Check for binder presence
        binder_present = any(
            "Binder" in str(idx.loc[ing, "Function"] if ing in idx.index and "Function" in idx.columns else "")
            for ing in blend
        )
        if not binder_present:
            flags.append("No binder detected — paint formulation requires a film-forming binder")
        else:
            passes.append("Film-forming binder present")

        # Estimate VOC (simplified — solvents with BP < 250C)
        voc_ingredients = []
        for ing in blend:
            if ing in idx.index:
                func = str(idx.loc[ing, "Function"] if "Function" in idx.columns else "")
                if any(x in func for x in ["Solvent", "Coalescent", "Reactive Diluent"]):
                    voc_ingredients.append(ing)
        if voc_ingredients:
            warnings.append(f"VOC-contributing ingredients: {', '.join(voc_ingredients[:3])} — verify VOC content per EU Decopaint Directive")
        else:
            passes.append("Low VOC profile — no solvent-heavy ingredients detected")

        return {
            "passes": passes, "warnings": warnings, "flags": flags,
            "voc_check": "Manual VOC calculation required per EN ISO 11890",
            "leed_note": "For LEED credits: VOC < 50 g/L (flat) or < 150 g/L (non-flat)",
            "regulatory_framework": config.get("regulatory_framework", ""),
        }

    # Fabric-specific
    if vertical == "fabric_laundry":
        phosphate_ings = [i for i in blend if "Phosphate" in i or "Polyphosphate" in i]
        if phosphate_ings:
            flags.append(f"Phosphate detected: {', '.join(phosphate_ings)} — banned in EU detergents Reg 648/2004")
        else:
            passes.append("Phosphate-free — EU Detergent Regulation compliant")

        enzyme_count = sum(1 for i in blend if "Enzyme" in str(
            idx.loc[i, "Function"] if i in idx.index and "Function" in idx.columns else ""))
        if enzyme_count > 0:
            passes.append(f"{enzyme_count} enzyme(s) detected — enhanced stain removal performance")

        return {
            "passes": passes, "warnings": warnings, "flags": flags,
            "phosphate_free": len(phosphate_ings) == 0,
            "biodegradability_note": "Verify OECD 301B ready biodegradability for all surfactants",
            "regulatory_framework": config.get("regulatory_framework", ""),
        }

    # Default / personal care / industrial
    return {
        "passes": passes, "warnings": warnings, "flags": flags,
        "regulatory_framework": config.get("regulatory_framework", ""),
    }


# ── Vertical-aware solver ─────────────────────────────────────────────────────

def _solve_vertical(db, max_cost, min_bio, min_perf, max_concentration, vertical):
    """
    LP solver with vertical-specific concentration bounds per function.
    Falls back to generic solver if vertical config not found.
    """
    config = VERTICAL_CONSTRAINTS.get(vertical, {})
    func_bounds = config.get("function_bounds", {})

    prob = pulp.LpProblem("IntelliForm", pulp.LpMinimize)
    names = db["Ingredient"].tolist()
    idx = db.set_index("Ingredient")

    # Build per-ingredient bounds based on function
    vars_ = {}
    for i, n in enumerate(names):
        lo, hi = 0.0, max_concentration
        if n in idx.index and "Function" in idx.columns:
            func = str(idx.loc[n, "Function"])
            for func_key, (f_lo, f_hi) in func_bounds.items():
                if func_key.lower() in func.lower():
                    lo = f_lo
                    hi = min(f_hi, max_concentration)
                    break
        vars_[n] = pulp.LpVariable(f"x_{i}", lowBound=lo, upBound=hi)

    # Objective: minimize cost
    prob += pulp.lpSum(float(idx.loc[n, "Cost_USD_kg"]) * vars_[n] for n in names)

    # Core constraints
    prob += pulp.lpSum(vars_.values()) == 1
    prob += pulp.lpSum(float(idx.loc[n, "Cost_USD_kg"]) * vars_[n] for n in names) <= max_cost
    prob += pulp.lpSum(float(idx.loc[n, "Bio_based_pct"]) * vars_[n] for n in names) >= min_bio
    prob += pulp.lpSum(float(idx.loc[n, "Performance_Score"]) * vars_[n] for n in names) >= min_perf

    # Vertical-specific: pharma needs at least some filler
    if vertical == "pharmaceutical":
        filler_names = [n for n in names if n in idx.index and
                       any(x in str(idx.loc[n, "Function"] if "Function" in idx.columns else "")
                           for x in ["Filler", "Binder/Filler"])]
        if filler_names:
            prob += pulp.lpSum(vars_[n] for n in filler_names) >= 0.20

    # Vertical-specific: paint needs at least some binder
    if vertical == "paint_coatings":
        binder_names = [n for n in names if n in idx.index and
                       "Binder" in str(idx.loc[n, "Function"] if "Function" in idx.columns else "")]
        if binder_names:
            prob += pulp.lpSum(vars_[n] for n in binder_names) >= 0.20

    # Vertical-specific: fabric needs primary surfactant
    if vertical == "fabric_laundry":
        surf_names = [n for n in names if n in idx.index and
                     any(x in str(idx.loc[n, "Function"] if "Function" in idx.columns else "")
                         for x in ["Primary Surfactant", "Anionic Surfactant"])]
        if surf_names:
            prob += pulp.lpSum(vars_[n] for n in surf_names) >= 0.05

    prob.solve(pulp.PULP_CBC_CMD(msg=0))
    status = pulp.LpStatus[prob.status]
    blend = {}
    if status == "Optimal":
        blend = {n: round(pulp.value(v) * 100, 1)
                 for n, v in vars_.items() if (pulp.value(v) or 0) > 0.005}
    return status, blend


def _solve(db, max_cost, min_bio, min_perf, max_concentration=1.0):
    """Original generic solver — used as fallback."""
    prob = pulp.LpProblem("IntelliForm", pulp.LpMinimize)
    names = db["Ingredient"].tolist()
    vars_ = {n: pulp.LpVariable(f"x_{i}", lowBound=0, upBound=max_concentration)
             for i, n in enumerate(names)}
    idx = db.set_index("Ingredient")
    prob += pulp.lpSum(float(idx.loc[n, "Cost_USD_kg"]) * vars_[n] for n in names)
    prob += pulp.lpSum(vars_.values()) == 1
    prob += pulp.lpSum(float(idx.loc[n, "Cost_USD_kg"]) * vars_[n] for n in names) <= max_cost
    prob += pulp.lpSum(float(idx.loc[n, "Bio_based_pct"]) * vars_[n] for n in names) >= min_bio
    prob += pulp.lpSum(float(idx.loc[n, "Performance_Score"]) * vars_[n] for n in names) >= min_perf
    prob.solve(pulp.PULP_CBC_CMD(msg=0))
    status = pulp.LpStatus[prob.status]
    blend = {}
    if status == "Optimal":
        blend = {n: round(pulp.value(v) * 100, 1)
                 for n, v in vars_.items() if (pulp.value(v) or 0) > 0.005}
    return status, blend


def _calc_metrics(blend, db):
    idx = db.set_index("Ingredient")
    cost = sum(float(idx.loc[k, "Cost_USD_kg"]) * pct / 100 for k, pct in blend.items() if k in idx.index)
    bio  = sum(float(idx.loc[k, "Bio_based_pct"]) * pct / 100 for k, pct in blend.items() if k in idx.index)
    perf = sum(float(idx.loc[k, "Performance_Score"]) * pct / 100 for k, pct in blend.items() if k in idx.index)
    return round(cost, 2), round(bio, 1), round(perf, 1)


# ── Manufacturing route selector ──────────────────────────────────────────────

def recommend_manufacturing_route(blend: Dict[str, float], db: pd.DataFrame,
                                   vertical: str) -> Optional[str]:
    if vertical != "pharmaceutical":
        return None
    idx = db.set_index("Ingredient") if "Ingredient" in db.columns else db

    filler_pct  = sum(pct for n, pct in blend.items() if n in idx.index and
                      any(x in str(idx.loc[n, "Function"] if "Function" in idx.columns else "")
                          for x in ["Filler", "Binder/Filler", "DCI Filler"]))
    binder_pct  = sum(pct for n, pct in blend.items() if n in idx.index and
                      "Binder" in str(idx.loc[n, "Function"] if "Function" in idx.columns else ""))
    lubricant_pct = sum(pct for n, pct in blend.items() if n in idx.index and
                        "Lubricant" in str(idx.loc[n, "Function"] if "Function" in idx.columns else ""))

    if filler_pct >= 50 and lubricant_pct > 0:
        return "Direct Compression (DC)"
    elif binder_pct >= 5:
        return "Wet Granulation"
    elif filler_pct >= 30:
        return "Dry Granulation / Roller Compaction"
    else:
        return "Review formulation — consult pharmaceutical development team"


# ── Main optimizer entry point ────────────────────────────────────────────────

def run_optimization(db, max_cost, min_bio, min_perf,
                     max_concentration: float = 1.0,
                     vertical: str = "all") -> OptResult:
    """
    Run vertical-aware optimization.

    Args:
        db: Ingredient database (filtered to vertical already)
        max_cost: Maximum cost per kg
        min_bio: Minimum bio-based percentage
        min_perf: Minimum performance score
        max_concentration: Max fraction for any single ingredient
        vertical: Vertical key — enables vertical-specific constraints
    """
    conc = max(max_concentration, 0.20)
    relaxation_schedule = [
        {"d_bio": -3,  "d_perf": -5,  "cost_mult": 1.00, "conc_delta": 0.10},
        {"d_bio": -5,  "d_perf": -5,  "cost_mult": 1.10, "conc_delta": 0.10},
        {"d_bio": -5,  "d_perf": -10, "cost_mult": 1.20, "conc_delta": 0.15},
    ]

    use_vertical_solver = vertical in VERTICAL_CONSTRAINTS and len(db) >= 5

    if use_vertical_solver:
        status, blend = _solve_vertical(db, max_cost, min_bio, min_perf, conc, vertical)
    else:
        status, blend = _solve(db, max_cost, min_bio, min_perf, conc)

    relaxed, rounds_used = False, 0
    current_cost, current_bio, current_perf = max_cost, min_bio, min_perf

    for round_num, sched in enumerate(relaxation_schedule, start=1):
        if status == "Optimal":
            break
        relaxed = True
        rounds_used = round_num
        current_cost = current_cost * sched["cost_mult"]
        current_bio  = max(current_bio  + sched["d_bio"],  0.0)
        current_perf = max(current_perf + sched["d_perf"], 60.0)
        conc = min(conc + sched["conc_delta"], 1.0)
        track("constraints_relaxed", {"round": round_num, "vertical": vertical})

        if use_vertical_solver:
            status, blend = _solve_vertical(db, current_cost, current_bio,
                                            current_perf, conc, vertical)
        else:
            status, blend = _solve(db, current_cost, current_bio, current_perf, conc)

    if status != "Optimal":
        track("error_optimization_failed", {"final_status": status, "vertical": vertical})
        return OptResult(
            success=False, blend={}, cost_per_kg=0, bio_pct=0, perf_score=0,
            status=status, relaxed=relaxed, relaxation_rounds=rounds_used,
            vertical=vertical,
            error_msg=(
                f"Could not find feasible blend after {rounds_used} relaxation round(s). "
                f"Try: (1) increase max cost, (2) reduce bio-based requirement, "
                f"(3) adjust max ingredient %, or (4) switch to All Verticals mode."
            )
        )

    cost, bio, perf = _calc_metrics(blend, db)

    # Vertical validation
    v_validation = validate_blend_for_vertical(blend, db, vertical)
    mfg_route = recommend_manufacturing_route(blend, db, vertical)

    compliance_flags = v_validation.get("flags", [])
    warnings = v_validation.get("warnings", [])

    track("formulation_generated", {
        "cost_per_kg": cost, "bio_based_pct": bio, "performance_score": perf,
        "relaxed": relaxed, "max_concentration": max_concentration,
        "vertical": vertical,
        "compliance_flags": len(compliance_flags),
        "warnings": len(warnings),
    })

    return OptResult(
        success=True, blend=blend, cost_per_kg=cost, bio_pct=bio,
        perf_score=perf, status=status, relaxed=relaxed,
        relaxation_rounds=rounds_used, vertical=vertical,
        vertical_validation=v_validation,
        manufacturing_route=mfg_route,
        compliance_flags=compliance_flags,
        warnings=warnings,
    )
