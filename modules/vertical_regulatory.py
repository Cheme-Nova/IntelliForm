"""
modules/vertical_regulatory.py
Vertical-specific regulatory intelligence for IntelliForm v1.1.

Each vertical has its own regulatory report:
  - pharmaceutical: ICH Q8/Q9/Q10, USP-NF, FDA 21 CFR, Ph.Eur
  - food: FDA GRAS, EU Food Additives, EFSA, Codex Alimentarius
  - agricultural: EPA FIFRA, EU 1107/2009, OMRI, REI/PHI
  - paint_coatings: EU Decopaint Directive, VOC Reg, LEED
  - fabric_laundry: EU Detergent Reg 648/2004, EPA DfE
  - industrial: REACH, EPA Safer Choice, NSF A1
  - personal_care: EU Cosmetics 1223/2009, COSMOS, REACH
"""
from dataclasses import dataclass, field
from typing import Dict, List, Optional
import pandas as pd


@dataclass
class VerticalRegulatoryReport:
    vertical: str
    overall_status: str          # ✅ Clear / ⚠️ Review Required / ❌ Blocked
    framework: str               # Primary regulatory framework
    passes: List[str]            # Items that passed
    warnings: List[str]          # Items needing review
    flags: List[str]             # Blocking issues
    certifications: List[str]    # Achievable certifications
    prohibited: List[str]        # Ingredients that violate regulations
    notes: List[str]             # General regulatory notes
    per_ingredient: List[dict]   # Per-ingredient regulatory detail


# ── Pharma regulatory lists ───────────────────────────────────────────────────

ICH_Q3A_SOLVENTS_CLASS1 = [
    "benzene", "carbon tetrachloride", "dichloroethane", "dichloromethane",
]
ICH_Q3A_SOLVENTS_CLASS2 = [
    "acetonitrile", "methanol", "n-hexane", "toluene", "chloroform",
]
ICH_SOLVENTS_CLASS3 = [
    "ethanol", "acetone", "isopropyl alcohol", "ethyl acetate", "propylene glycol",
    "water", "glycerol",
]

PHARMA_BANNED = [
    "benzalkonium chloride",   # restricted in some routes
    "thiomersal",              # mercury-based, restricted
    "propylene glycol",        # limited in pediatric oral
]

USP_EXCIPIENTS = [
    "microcrystalline cellulose", "lactose monohydrate", "magnesium stearate",
    "hydroxypropyl methylcellulose", "povidone", "croscarmellose sodium",
    "talc", "silicon dioxide", "stearic acid", "sodium starch glycolate",
    "gelatin", "titanium dioxide", "calcium phosphate", "corn starch",
    "mannitol", "sorbitol", "sucrose", "glycerin", "propylene glycol",
    "ethanol", "water for injection",
]

# ── Food regulatory lists ─────────────────────────────────────────────────────

EU_FOOD_PROHIBITED = [
    "titanium dioxide",   # EU banned food additive E171 since 2022
]
MAJOR_ALLERGENS = [
    "wheat", "gluten", "milk", "casein", "whey", "egg", "albumen",
    "peanut", "soy", "tree nut", "fish", "shellfish", "sesame",
    "mustard", "celery", "lupin", "molluscs", "sulphites", "sulfites",
]
CLEAN_LABEL_AVOID = [
    "sodium nitrite", "sulfur dioxide", "sodium erythorbate", "bht",
    "bha", "tbhq", "sodium benzoate", "potassium sorbate",
    "carrageenan", "titanium dioxide", "high fructose",
]
ORGANIC_ALLOWED = [
    "lecithin", "pectin", "agar", "carob bean gum", "guar gum",
    "xanthan gum", "citric acid", "ascorbic acid", "tocopherols",
    "rosemary extract", "potassium sorbate",
]

# ── Agricultural regulatory lists ─────────────────────────────────────────────

OMRI_LISTED = [
    "neem oil", "pyrethrin", "spinosad", "bacillus thuringiensis",
    "beauveria bassiana", "kaolin", "copper hydroxide", "sulfur",
    "hydrogen peroxide", "citric acid", "ethanol", "potassium bicarbonate",
    "insecticidal soap", "diatomaceous earth", "iron phosphate",
    "bacillus subtilis", "trichoderma",
]
EU_BANNED_AGRI = [
    "chlorpyrifos", "paraquat", "neonicotinoid",
]

# ── Paint regulatory lists ────────────────────────────────────────────────────

SVHC_PIGMENTS = [
    "lead chromate", "cadmium", "hexavalent chromium",
    "cobalt dryer",  # cobalt compounds SVHC
]
HIGH_VOC_SOLVENTS = [
    "toluene", "xylene", "benzene", "styrene", "ethylbenzene",
    "methyl ethyl ketone", "n-butyl acetate",
]

# ── Fabric regulatory lists ───────────────────────────────────────────────────

EU_BANNED_FABRIC = [
    "alkylphenol ethoxylate", "nonylphenol", "tributyltin",
    "perfluorooctane", "pfas",
]
PHOSPHATE_RESTRICTED = [
    "sodium tripolyphosphate", "trisodium phosphate",
    "sodium hexametaphosphate", "tetrasodium pyrophosphate",
]


# ── Core regulatory engine ────────────────────────────────────────────────────

def generate_vertical_regulatory_report(
    blend: Dict[str, float],
    db: pd.DataFrame,
    vertical: str,
) -> VerticalRegulatoryReport:
    """Generate a vertical-specific regulatory report for a blend."""

    idx = db.set_index("Ingredient") if "Ingredient" in db.columns else db
    passes, warnings, flags, notes = [], [], [], []
    certifications = []
    prohibited = []
    per_ingredient = []

    ing_names_lower = {ing.lower(): ing for ing in blend}

    if vertical == "pharmaceutical":
        framework = "ICH Q8/Q9/Q10 · USP-NF · FDA 21 CFR 211 · Ph.Eur"
        notes.append("All excipients should comply with USP-NF or Ph.Eur monographs")
        notes.append("ICH Q3C solvent residue limits apply to manufacturing solvents")
        notes.append("Excipient qualification per ICH Q7 required for API contact materials")

        for ing, pct in blend.items():
            ing_l = ing.lower()
            status = "✅"
            reg_notes = []

            # Check against USP list
            if any(u in ing_l for u in USP_EXCIPIENTS):
                reg_notes.append("USP/NF listed excipient")
                passes.append(f"{ing}: USP/NF compendial grade available")
            else:
                reg_notes.append("Verify compendial status")
                warnings.append(f"{ing}: Confirm USP/NF or Ph.Eur monograph exists")
                status = "⚠️"

            # ICH Q3C solvent check
            if any(s in ing_l for s in ICH_Q3A_SOLVENTS_CLASS1):
                flags.append(f"{ing}: ICH Q3C Class 1 solvent — should be avoided")
                prohibited.append(ing)
                status = "❌"
            elif any(s in ing_l for s in ICH_Q3A_SOLVENTS_CLASS2):
                warnings.append(f"{ing}: ICH Q3C Class 2 solvent — strict limits apply")
                status = "⚠️"

            per_ingredient.append({
                "ingredient": ing, "pct": pct, "status": status,
                "notes": "; ".join(reg_notes) if reg_notes else "Standard excipient",
                "ich_class": "Class 3" if any(s in ing_l for s in ICH_SOLVENTS_CLASS3) else "—",
            })

        if not flags:
            passes.append("No ICH Class 1 solvents detected")
        certifications = ["GMP Manufacturing Ready", "ICH Q8 Compliant Formulation Design"]
        overall = "✅ Clear" if not flags else ("❌ Blocked" if prohibited else "⚠️ Review Required")

    elif vertical == "food":
        framework = "FDA GRAS · EU Food Additives Reg 1333/2008 · EFSA · Codex Alimentarius"
        notes.append("All additives must have GRAS status (FDA) or approved EU food additive number")
        notes.append("Allergen labeling mandatory per EU 1169/2011 and FDA FALCPA")
        notes.append("Organic certification requires USDA NOP / EU 2018/848 approved ingredients")

        allergens_found = []
        clean_label_issues = []
        eu_banned_found = []

        for ing, pct in blend.items():
            ing_l = ing.lower()
            status = "✅"
            reg_notes = []

            # EU banned check
            if any(b in ing_l for b in EU_FOOD_PROHIBITED):
                flags.append(f"{ing}: BANNED in EU food applications since 2022")
                eu_banned_found.append(ing)
                prohibited.append(ing)
                status = "❌"

            # Allergen check
            allergen_match = [a for a in MAJOR_ALLERGENS if a in ing_l]
            if allergen_match:
                allergens_found.append(ing)
                reg_notes.append(f"ALLERGEN: {allergen_match[0].upper()}")
                warnings.append(f"{ing}: Allergen declaration required ({allergen_match[0]})")
                status = "⚠️"

            # Clean label check
            if any(c in ing_l for c in CLEAN_LABEL_AVOID):
                clean_label_issues.append(ing)
                reg_notes.append("Not clean-label")

            # Organic check
            if any(o in ing_l for o in ORGANIC_ALLOWED):
                reg_notes.append("USDA Organic / EU Organic allowed")
            else:
                reg_notes.append("Verify organic certification status")

            per_ingredient.append({
                "ingredient": ing, "pct": pct, "status": status,
                "notes": "; ".join(reg_notes) if reg_notes else "Standard food additive",
                "gras": "Assumed — verify",
            })

        if allergens_found:
            flags.append(f"MANDATORY ALLERGEN LABELING: {', '.join(allergens_found)}")
        else:
            passes.append("No major allergens detected — allergen-free claim supportable")

        if not clean_label_issues:
            passes.append("Clean label compliant")
            certifications.append("Clean Label Claim Supportable")
        else:
            warnings.append(f"Non-clean-label ingredients: {', '.join(clean_label_issues)}")

        if not eu_banned_found:
            passes.append("No EU-banned food additives detected")
            certifications.append("EU Food Additive Compliant")

        certifications.append("GRAS Presumed — formal FDA GRAS notification recommended for novel uses")
        overall = "✅ Clear" if not flags else ("❌ Blocked" if prohibited else "⚠️ Review Required")

    elif vertical == "agricultural":
        framework = "EPA FIFRA · EU Reg 1107/2009 · OMRI Listed · NOP (USDA Organic)"
        notes.append("EPA registration required before commercial sale in the US")
        notes.append("EU Basic Substance or Low Risk Active may apply for bio-based ingredients")
        notes.append("Verify REI and PHI for each active ingredient on crop labels")
        notes.append("Tank mix compatibility testing required before commercial use")

        omri_listed_found = []
        eu_concern = []

        for ing, pct in blend.items():
            ing_l = ing.lower()
            status = "✅"
            reg_notes = []

            # OMRI check
            if any(o in ing_l for o in OMRI_LISTED):
                omri_listed_found.append(ing)
                reg_notes.append("OMRI Listed — organic approved")
                passes.append(f"{ing}: OMRI Listed")
            else:
                reg_notes.append("Verify OMRI/organic status")

            # EU banned check
            if any(b in ing_l for b in EU_BANNED_AGRI):
                flags.append(f"{ing}: Restricted or banned in EU agriculture")
                eu_concern.append(ing)
                status = "❌"

            # Potency warning
            if any(x in ing_l for x in ["plant growth regulator", "gibberellic", "cytokinin", "auxin"]):
                warnings.append(f"{ing}: PGR — EPA registration required, strict dosing required")
                status = "⚠️"

            per_ingredient.append({
                "ingredient": ing, "pct": pct, "status": status,
                "notes": "; ".join(reg_notes) if reg_notes else "Verify EPA/EU registration",
                "omri": "✅" if any(o in ing_l for o in OMRI_LISTED) else "Verify",
                "rei": "See label", "phi": "See label",
            })

        if omri_listed_found:
            certifications.append(f"OMRI Organic Compliant — {len(omri_listed_found)} OMRI ingredients")
        if not eu_concern:
            passes.append("No EU-banned active substances detected")

        certifications.append("EPA FIFRA — Registration required before commercial sale")
        certifications.append("EU Basic Substance review may apply for food-grade actives")
        overall = "✅ Clear" if not flags else ("❌ Blocked" if prohibited else "⚠️ Review Required")

    elif vertical == "paint_coatings":
        framework = "EU Decopaint Directive 2004/42/EC · VOC Reg EC 2004/42 · REACH · LEED v4"
        notes.append("VOC content must comply with EU Decopaint Directive phase 2 limits")
        notes.append("SVHC ingredients require REACH Article 33 disclosure >0.1% w/w")
        notes.append("LEED v4 credit requires VOC <50 g/L for architectural flat paints")
        notes.append("EN ISO 11890-1/2 for VOC determination")

        svhc_found = []
        high_voc_found = []
        binder_count = 0

        for ing, pct in blend.items():
            ing_l = ing.lower()
            status = "✅"
            reg_notes = []

            # SVHC check
            if any(s in ing_l for s in SVHC_PIGMENTS):
                flags.append(f"{ing}: SVHC / restricted substance — REACH disclosure required")
                svhc_found.append(ing)
                prohibited.append(ing)
                status = "❌"

            # VOC check
            if any(v in ing_l for v in HIGH_VOC_SOLVENTS):
                high_voc_found.append(ing)
                warnings.append(f"{ing}: High-VOC solvent — check EU Decopaint limit for product category")
                status = "⚠️"

            if "binder" in ing_l or "resin" in ing_l or "alkyd" in ing_l:
                binder_count += 1
                reg_notes.append("Film former — check SVHC registration")

            per_ingredient.append({
                "ingredient": ing, "pct": pct, "status": status,
                "notes": "; ".join(reg_notes) if reg_notes else "Standard coatings ingredient",
                "reach": "Verify ECHA registration",
                "voc_contributor": "Yes" if any(v in ing_l for v in HIGH_VOC_SOLVENTS) else "Low",
            })

        if not svhc_found:
            passes.append("No SVHC substances detected")
            certifications.append("REACH Compliant — no restricted substances")
        if not high_voc_found:
            passes.append("Low-VOC profile — likely EU Decopaint compliant")
            certifications.append("Low-VOC / Zero-VOC Claim Supportable")
            certifications.append("LEED v4 Low-Emitting Materials Credit")
        if binder_count > 0:
            passes.append(f"Film-forming system present ({binder_count} binder(s))")

        overall = "✅ Clear" if not flags else ("❌ Blocked" if prohibited else "⚠️ Review Required")

    elif vertical == "fabric_laundry":
        framework = "EU Detergent Regulation 648/2004 · REACH · EPA DfE · Cradle to Cradle"
        notes.append("All surfactants must meet OECD 301B ready biodegradability (>60% in 28 days)")
        notes.append("Fragrances must comply with IFRA guidelines and EU Annex III")
        notes.append("Phosphate-free mandatory for EU household detergents since 2013")
        notes.append("Ingredient disclosure required on product label and online")

        phosphate_found = []
        restricted_found = []
        enzyme_count = 0

        for ing, pct in blend.items():
            ing_l = ing.lower()
            status = "✅"
            reg_notes = []

            # Phosphate check
            if any(p in ing_l for p in PHOSPHATE_RESTRICTED):
                flags.append(f"{ing}: PHOSPHATE — banned in EU household detergents since 2013")
                phosphate_found.append(ing)
                prohibited.append(ing)
                status = "❌"

            # EU banned check
            if any(b in ing_l for b in EU_BANNED_FABRIC):
                flags.append(f"{ing}: Restricted substance in fabric care")
                restricted_found.append(ing)
                status = "❌"

            # Enzyme detection
            if any(x in ing_l for x in ["protease", "lipase", "amylase", "cellulase",
                                         "savinase", "natalase", "lipolase", "carezyme",
                                         "everlase", "esperase"]):
                enzyme_count += 1
                reg_notes.append("Enzyme — safe handling practices required (respiratory sensitizer)")
                warnings.append(f"{ing}: Enzyme — occupational exposure limits apply during manufacturing")
                status = "⚠️" if status == "✅" else status

            per_ingredient.append({
                "ingredient": ing, "pct": pct, "status": status,
                "notes": "; ".join(reg_notes) if reg_notes else "Standard detergent ingredient",
                "biodeg": "Verify OECD 301B",
                "reach": "Verify ECHA",
            })

        if not phosphate_found:
            passes.append("Phosphate-free — EU Detergent Regulation compliant")
            certifications.append("EU Detergent Regulation 648/2004 Compliant")
        if enzyme_count > 0:
            certifications.append(f"Enzymatic Cleaning — {enzyme_count} enzyme(s) for enhanced stain removal")
        if not restricted_found:
            passes.append("No restricted substances detected")
            certifications.append("EPA Design for Environment (DfE) Eligible")

        overall = "✅ Clear" if not flags else ("❌ Blocked" if prohibited else "⚠️ Review Required")

    elif vertical == "industrial":
        framework = "REACH · EPA Safer Choice · NSF/ANSI 61 · NSF A1 (food-safe)"
        notes.append("SDS preparation required for all industrial formulations (GHS/CLP)")
        notes.append("EPA Safer Choice certification improves market access in institutional segment")
        notes.append("NSF A1 approval required for food-processing facility use")

        svhc_found = []

        for ing, pct in blend.items():
            ing_l = ing.lower()
            status = "✅"
            reg_notes = []

            if any(x in ing_l for x in ["benzene", "formaldehyde", "chromate", "asbestos"]):
                flags.append(f"{ing}: Restricted under REACH SVHC list")
                svhc_found.append(ing)
                prohibited.append(ing)
                status = "❌"

            if any(x in ing_l for x in ["isothiazolinone", "glutaraldehyde", "quaternary ammonium"]):
                warnings.append(f"{ing}: Biocide — EU BPR registration required for biocidal claims")
                status = "⚠️"

            if "epa safer" in ing_l or "bio-" in ing_l or "bio " in ing_l:
                reg_notes.append("EPA Safer Choice ingredient candidate")
                passes.append(f"{ing}: EPA Safer Choice eligible ingredient")

            per_ingredient.append({
                "ingredient": ing, "pct": pct, "status": status,
                "notes": "; ".join(reg_notes) if reg_notes else "Standard industrial ingredient",
                "reach": "Verify ECHA registration",
                "sds_required": "Yes",
            })

        if not svhc_found:
            passes.append("No SVHC substances detected")
            certifications.append("REACH Compliant")
            certifications.append("EPA Safer Choice Eligible — certification application recommended")

        overall = "✅ Clear" if not flags else ("❌ Blocked" if prohibited else "⚠️ Review Required")

    else:
        # Personal care / default
        framework = "EU Cosmetics Regulation 1223/2009 · COSMOS-standard · REACH"
        notes.append("EU Cosmetics Regulation requires safety assessment by qualified toxicologist")
        notes.append("INCI naming required on product labels")
        notes.append("Prohibited/restricted substances listed in Annexes II-VI of EU Reg 1223/2009")
        notes.append("Nanomaterials require special labeling and notification")

        restricted_found = []
        eu_prohibited = ["formaldehyde", "mercury", "lead", "arsenic", "hydroquinone"]

        for ing, pct in blend.items():
            ing_l = ing.lower()
            status = "✅"
            reg_notes = []

            if any(r in ing_l for r in eu_prohibited):
                flags.append(f"{ing}: Prohibited under EU Cosmetics Reg Annex II")
                restricted_found.append(ing)
                prohibited.append(ing)
                status = "❌"

            if "fragrance" in ing_l or "essential oil" in ing_l:
                warnings.append(f"{ing}: Fragrance allergen labeling required if >0.001% (EU 2023/1490)")
                status = "⚠️"

            if any(x in ing_l for x in ["titanium dioxide", "zinc oxide"]) and pct < 25:
                warnings.append(f"{ing}: Nanomaterial form requires EU notification and [nano] label")

            per_ingredient.append({
                "ingredient": ing, "pct": pct, "status": status,
                "notes": "; ".join(reg_notes) if reg_notes else "Verify EU Cosmetics Reg compliance",
                "inci_name": "Verify INCI name",
                "cosmos": "Check COSMOS approved list",
            })

        if not restricted_found:
            passes.append("No EU Annex II prohibited substances detected")
            certifications.append("EU Cosmetics Regulation 1223/2009 Compliant")

        cosmos_eligible = all(
            any(x in ing.lower() for x in ["bio", "natural", "plant", "extract", "oil"])
            or ing.lower() in ["water", "citric acid", "lactic acid", "sodium chloride"]
            for ing in blend
        )
        if cosmos_eligible:
            certifications.append("COSMOS Natural/Organic — eligibility check recommended")
        certifications.append("REACH — verify ECHA registration for all ingredients")

        overall = "✅ Clear" if not flags else ("❌ Blocked" if prohibited else "⚠️ Review Required")

    return VerticalRegulatoryReport(
        vertical=vertical,
        overall_status=overall,
        framework=framework,
        passes=passes,
        warnings=warnings,
        flags=flags,
        certifications=certifications,
        prohibited=prohibited,
        notes=notes,
        per_ingredient=per_ingredient,
    )
