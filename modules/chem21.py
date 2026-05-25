"""
modules/chem21.py
CHEM21 Solvent Greenness Scoring for IntelliForm.

Implements the pharmaceutical industry's standard solvent selection guide,
based on the CHEM21 project (Prat et al., Green Chemistry 2016, DOI 10.1039/C5GC01015G)
augmented with GlaxoSmithKline, Pfizer, and Sanofi solvent guides.

Tiers:
  1 — Recommended   (green)   : low hazard, low environmental impact
  2 — Problematic   (yellow)  : usable with care, substitution encouraged
  3 — Hazardous     (amber)   : significant concerns, substitution required where feasible
  4 — Highly Hazardous (red)  : restricted, phase-out recommended

Solvent Greenness Score (0-100): 100 = all Recommended, 0 = all Highly Hazardous
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

# ── CHEM21 Solvent Database ────────────────────────────────────────────────────
# Tier 1=Recommended (100pt), 2=Problematic (60pt), 3=Hazardous (25pt), 4=Highly Hazardous (0pt)
# Sources: CHEM21 (2016), GSK Solvent Guide (2016), Pfizer Green Chemistry (2018), Sanofi (2013)

_CHEM21_DB: Dict[str, Tuple[int, str]] = {
    # ── TIER 1 — RECOMMENDED ──────────────────────────────────────────────────
    "water":                   (1, "Ideal solvent. No VOC, non-toxic, non-flammable."),
    "ethanol":                 (1, "Renewable, low toxicity, widely available."),
    "ethanol (bio)":           (1, "Bio-based ethanol — preferred over petrochemical."),
    "isopropanol":             (1, "Low toxicity; reproductive concern at high exposure."),
    "isopropyl alcohol":       (1, "Low toxicity; reproductive concern at high exposure."),
    "1-propanol":              (1, "Low hazard profile; preferred over higher alcohols."),
    "n-propanol":              (1, "Low hazard profile."),
    "1-butanol":               (1, "Recommended; irritant but low systemic toxicity."),
    "n-butanol":               (1, "Recommended; irritant but low systemic toxicity."),
    "ethyl acetate":           (1, "Low toxicity, low boiling point, biodegradable."),
    "isopropyl acetate":       (1, "Preferred ester solvent; good biodegradation."),
    "acetone":                 (1, "Low toxicity, widely available, biodegradable."),
    "methyl ethyl ketone":     (1, "MEK — low acute toxicity, biodegradable."),
    "mek":                     (1, "Methyl ethyl ketone — low acute toxicity."),
    "heptane":                 (1, "Preferred alkane over hexane — lower neurotoxicity."),
    "n-heptane":               (1, "Preferred alkane — lower neurotoxicity than hexane."),
    "cyclohexane":             (1, "Moderate hazard — preferred over benzene/toluene."),
    "methylcyclohexane":       (1, "Moderate hazard, lower toxicity than cyclohexane."),
    "glycerol":                (1, "Bio-derived, non-toxic, biodegradable."),
    "glycerol (fabric)":       (1, "Bio-derived, non-toxic, biodegradable."),
    "propylene glycol":        (1, "Low toxicity, biodegradable, GRAS status."),
    "1,3-propanediol":         (1, "Bio-based (DuPont Zemea), preferred over PG."),
    "diethyl carbonate":       (1, "Preferred carbonates — low toxicity."),
    "dimethyl carbonate":      (1, "Low toxicity alternative to DCM/DMF."),
    "2-methylthf":             (1, "Bio-based from furfural, lower aquatic toxicity than THF."),
    "2-methyl tetrahydrofuran":(1, "Bio-based from furfural — recommended THF replacement."),
    "cyrene":                  (1, "Bio-derived from cellulose, highly recommended."),
    "glycofurol":              (1, "Low toxicity pharmaceutical grade solvent."),
    "d-limonene":              (1, "Bio-based terpene, low systemic toxicity."),
    "ethyl lactate":           (1, "Bio-derived, biodegradable, low toxicity."),
    "lactic acid":             (1, "Bio-derived acid, GRAS, biodegradable."),
    "levulinic acid":          (1, "Bio-based platform chemical, low toxicity."),
    "citric acid":             (1, "Bio-derived, GRAS, biodegradable."),
    "acetic acid":             (1, "Widely available, biodegradable; strong acid — handle with care."),
    "methanol":                (1, "Low boiling, low cost; flammable — keep away from ignition."),
    "tert-butanol":            (1, "Recommended cosolvent; mild CNS effects at high dose."),
    "t-butanol":               (1, "Recommended cosolvent."),

    # ── TIER 2 — PROBLEMATIC ──────────────────────────────────────────────────
    "dimethyl sulfoxide":      (2, "DMSO: penetration enhancer — may carry contaminants through skin."),
    "dmso":                    (2, "Penetration enhancer — may carry contaminants through skin."),
    "acetonitrile":            (2, "Moderate toxicity, high VOC; substitute where feasible."),
    "methyl tert-butyl ether": (2, "MTBE: groundwater contamination risk; use 2-MeTHF instead."),
    "mtbe":                    (2, "Groundwater contamination risk; use 2-MeTHF."),
    "diethyl ether":           (2, "Extreme flammability, peroxide formation risk."),
    "tetrahydrofuran":         (2, "THF: moderate toxicity, peroxide formation; use 2-MeTHF."),
    "thf":                     (2, "Peroxide formation risk; use 2-methylTHF instead."),
    "xylene":                  (2, "Moderate toxicity, high VOC; substitute preferred."),
    "xylenes":                 (2, "Moderate toxicity, high VOC."),
    "toluene":                 (2, "CNS effects, reproductive concern — minimize use."),
    "dichloromethane":         (2, "DCM: CARCINOGEN (IARC 2A); substitute required in pharma post-2025."),
    "dcm":                     (2, "CARCINOGEN (IARC 2A); use Cyrene or EtOAc instead."),
    "methylene chloride":      (2, "CARCINOGEN (IARC 2A); phase-out mandated by EU."),
    "ethylene glycol":         (2, "Low acute toxicity but renal toxicity at high dose."),
    "dimethylformamide":       (2, "DMF: reproductive toxin (Category 1B) — restricted."),
    "dmf":                     (2, "Reproductive toxin — use DMAc alternatives."),
    "dimethylacetamide":       (2, "DMAc: reproductive toxin — minimize, substitute NMP alternatives."),
    "diethylene glycol":       (2, "Metabolizes to glycolic acid; renal toxin — substitution preferred."),
    "butyl acetate":           (2, "Moderate VOC; biodegradable but irritant at high concentration."),
    "isobutyl acetate":        (2, "Moderate VOC; biodegradable but irritant."),
    "isophorone":              (2, "Skin/eye irritant; moderate aquatic toxicity."),
    "propylene carbonate":     (2, "Moderate toxicity; degradable — use dimethyl carbonate where possible."),
    "2-butanol":               (2, "Mild CNS effects; substitution preferred with 1-butanol."),
    "isobutanol":              (2, "CNS effects at high dose; fire hazard."),
    "2-butoxyethanol":         (2, "Blood toxicity concern; restricted in many cosmetic applications."),
    "polyethylene glycol":     (2, "PEG: oligomeric contamination risks; dermal concerns at small MW."),
    "peg-400":                 (2, "Oligomeric contamination risk; COSMOS restricts some grades."),

    # ── TIER 3 — HAZARDOUS ────────────────────────────────────────────────────
    "n-methyl-2-pyrrolidone":  (3, "NMP: reproductive toxin (CMR Cat 1B); REACH restricted."),
    "nmp":                     (3, "Reproductive toxin; REACH authorisation required."),
    "chloroform":              (3, "Hepatotoxin, suspected carcinogen (IARC 2B) — prohibited in cosmetics."),
    "n-hexane":                (3, "Neurotoxin (hexane neuropathy via 2,5-hexanedione); substitute with heptane."),
    "hexane":                  (3, "Neurotoxin — replace with heptane in all applications."),
    "benzene":                 (3, "CARCINOGEN (IARC 1); no legitimate use as solvent."),
    "pyridine":                (3, "Reproductive toxin, high acute toxicity; substitution required."),
    "dioxane":                 (3, "1,4-Dioxane: CARCINOGEN (IARC 2B); contaminant in PEG ethoxylation."),
    "1,4-dioxane":             (3, "CARCINOGEN; restricted in cosmetics products as trace impurity."),
    "styrene":                 (3, "Suspected carcinogen; VOC; neurotoxin at chronic exposure."),
    "dioxolane":               (3, "Moderate to high toxicity; substitution recommended."),
    "dibromomethane":          (3, "Halogenated solvent; substitution required."),
    "nitromethane":            (3, "Explosive risk at high concentrations; highly toxic."),
    "diethanolamine":          (3, "DEA: contamination with nitrosamines; COSMOS restricted."),
    "triethanolamine":         (3, "TEA: nitrosamine contamination risk; pH-adjust alternatives preferred."),

    # ── TIER 4 — HIGHLY HAZARDOUS ────────────────────────────────────────────
    "carbon tetrachloride":    (4, "CARCINOGEN (IARC 1); hepatotoxin; ozone depleting substance."),
    "chlorobenzene":           (4, "High toxicity, environmental persistence — prohibited."),
    "carbon disulfide":        (4, "Highly toxic; neurotoxin; reproductive toxin."),
    "formaldehyde":            (4, "CARCINOGEN (IARC 1); restricted to trace levels in all sectors."),
    "nitrobenzene":            (4, "Highly toxic; methemoglobin-forming — prohibited as solvent."),
    "aniline":                 (4, "Methemoglobin-forming agent; carcinogen; prohibited."),
    "diisocyanates":           (4, "Sensitisers; severe respiratory hazard — restricted."),
    "trichloroethylene":       (4, "CARCINOGEN (IARC 1A); EU REACH restricted entry in cosmetics."),
    "perchloroethylene":       (4, "CARCINOGEN — REACH restricted; phase-out mandated."),
    "methyl iodide":           (4, "Carcinogen, alkylating agent — no cosmetic/food use."),
    "epichlorohydrin":         (4, "Carcinogen, mutagen, reproductive toxin."),
    "ethylene oxide":          (4, "CARCINOGEN (IARC 1) — restricted as gas sterilant only."),
    "propylene oxide":         (4, "CARCINOGEN (IARC 2B); mutagen."),
    "allyl chloride":          (4, "Neurotoxin, high reactivity — no legitimate solvent use."),
}

# Keyword matching for fuzzy ingredient name lookup
_SOLVENT_KEYWORDS = {
    "ethanol": "ethanol",
    "bio-ethanol": "ethanol (bio)",
    "isopropanol": "isopropanol",
    "ipa ": "isopropanol",
    "acetone": "acetone",
    "heptane": "heptane",
    "cyclohexane": "cyclohexane",
    "toluene": "toluene",
    "xylene": "xylene",
    "dmso": "dimethyl sulfoxide",
    "thf": "tetrahydrofuran",
    "dcm": "dichloromethane",
    "nmp": "n-methyl-2-pyrrolidone",
    "dmf": "dimethylformamide",
    "hexane": "hexane",
    "n-hexane": "n-hexane",
    "d-limonene": "d-limonene",
    "glycerol": "glycerol",
    "propylene glycol": "propylene glycol",
}

# Tier score mapping
_TIER_SCORES = {1: 100, 2: 60, 3: 25, 4: 0}
_TIER_LABELS = {1: "Recommended", 2: "Problematic", 3: "Hazardous", 4: "Highly Hazardous"}
_TIER_COLORS = {1: "green", 2: "yellow", 3: "amber", 4: "red"}


@dataclass
class SolventGreenness:
    """CHEM21 solvent assessment for a single ingredient."""
    name: str
    tier: int                    # 1–4
    tier_label: str              # "Recommended" etc.
    score: float                 # 0–100
    color: str                   # green/yellow/amber/red
    rationale: str
    substitute: Optional[str] = None


@dataclass
class BlendSolventProfile:
    """Aggregated CHEM21 profile for a full blend."""
    weighted_score: float          # 0–100, normalized to assessed solvents only
    grade: str                     # A+ / A / B / C / D
    coverage_pct: float            # % of blend weight covered by CHEM21 database
    solvents_assessed: List[SolventGreenness]
    worst_tier: int                # 1–4
    alerts: List[str]              # any tier-3/4 warnings
    substitution_suggestions: List[str]


# Non-solvent ingredient patterns to skip (food actives, functional ingredients, etc.)
_NON_SOLVENT_PATTERNS = {
    "lecithin", "sucrose", "fructose", "starch", "cellulose", "pectin",
    "gum", "enzyme", "protein", "peptide", "extract", "powder", "flour",
    "salt", "chloride", "sulfate", "phosphate", "carbonate", "bicarbonate",
    "hydroxide", "oxide", "silica", "titanium", "zinc", "iron", "calcium",
    "pigment", "wax", "ester (", "paraffin", "fragrance", "dye",
    "preservative", "benzoate", "sorbate", "propionate",
}


def _lookup_solvent(name: str) -> Optional[Tuple[int, str]]:
    """Look up a solvent in the CHEM21 database by ingredient name."""
    key = name.lower().strip()
    # Skip obvious non-solvents
    if any(pat in key for pat in _NON_SOLVENT_PATTERNS):
        return None
    if key in _CHEM21_DB:
        return _CHEM21_DB[key]
    # Fuzzy keyword match
    for keyword, canonical in _SOLVENT_KEYWORDS.items():
        if keyword in key:
            return _CHEM21_DB.get(canonical)
    # Partial match — only if the ingredient name starts with a known key
    # (avoids false matches like "citric acid" → "acetic acid")
    for db_key, val in _CHEM21_DB.items():
        if key.startswith(db_key) or db_key.startswith(key):
            return val
    return None


def assess_solvent(name: str) -> Optional[SolventGreenness]:
    """Return CHEM21 greenness assessment for a single ingredient name."""
    result = _lookup_solvent(name)
    if result is None:
        return None
    tier, rationale = result
    return SolventGreenness(
        name=name,
        tier=tier,
        tier_label=_TIER_LABELS[tier],
        score=float(_TIER_SCORES[tier]),
        color=_TIER_COLORS[tier],
        rationale=rationale,
    )


def score_blend_solvents(blend: dict) -> BlendSolventProfile:
    """
    Compute CHEM21 solvent greenness profile for a full blend.

    Args:
        blend: {ingredient_name: weight_percent}

    Returns:
        BlendSolventProfile with weighted score and alerts.
    """
    assessed = []
    total_solvent_pct = 0.0
    weighted_score = 0.0
    worst_tier = 1
    alerts = []
    suggestions = []

    for ingredient, pct in blend.items():
        sg = assess_solvent(ingredient)
        if sg is None:
            continue
        assessed.append(sg)
        total_solvent_pct += pct
        weighted_score += sg.score * (pct / 100.0)
        if sg.tier > worst_tier:
            worst_tier = sg.tier
        if sg.tier >= 3:
            alerts.append(f"{ingredient} ({sg.tier_label}): {sg.rationale[:80]}")
        if sg.tier >= 2:
            _add_suggestion(suggestions, sg)

    if not assessed:
        return BlendSolventProfile(
            weighted_score=75.0,
            grade="B",
            coverage_pct=0.0,
            solvents_assessed=[],
            worst_tier=1,
            alerts=[],
            substitution_suggestions=[],
        )

    # Normalize score to ONLY the assessed solvent fraction.
    # e.g. 100% Tier-1 ethanol at 25% of blend → score = 100, not 25.
    final_score = round((weighted_score / (total_solvent_pct / 100.0)), 1) if total_solvent_pct > 0 else 75.0

    return BlendSolventProfile(
        weighted_score=final_score,
        grade=_grade_solvent(final_score),
        coverage_pct=round(total_solvent_pct, 1),
        solvents_assessed=assessed,
        worst_tier=worst_tier,
        alerts=alerts,
        substitution_suggestions=list(set(suggestions)),
    )


def _add_suggestion(suggestions: list, sg: SolventGreenness):
    if "hexane" in sg.name.lower() and "heptane" not in sg.name.lower():
        suggestions.append("Replace n-hexane with heptane (lower neurotoxicity, CHEM21 Tier 1)")
    elif "dcm" in sg.name.lower() or "dichloromethane" in sg.name.lower():
        suggestions.append("Replace DCM with Cyrene or ethyl acetate (CHEM21 Recommended)")
    elif "dmf" in sg.name.lower():
        suggestions.append("Replace DMF with dimethyl carbonate or Cyrene (CHEM21 Tier 1)")
    elif "nmp" in sg.name.lower():
        suggestions.append("Replace NMP with Cyrene, GVL, or DMSO (lower CMR concern)")
    elif "toluene" in sg.name.lower():
        suggestions.append("Replace toluene with ethyl acetate or methylcyclohexane (CHEM21 Tier 1)")
    elif "thf" in sg.name.lower():
        suggestions.append("Replace THF with 2-methylTHF (bio-based, lower aquatic toxicity)")


def _grade_solvent(score: float) -> str:
    if score >= 90: return "A+"
    if score >= 75: return "A"
    if score >= 55: return "B"
    if score >= 30: return "C"
    return "D"
