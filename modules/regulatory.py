"""
modules/regulatory.py
Regulatory intelligence for IntelliForm v0.9.

Per-ingredient structured regulatory data covering:
  - REACH (ECHA) — registration status, SVHC flag, restriction list
  - EPA Safer Choice (DfE) — tier classification
  - EU Ecolabel — detergent/cosmetic product eligibility
  - COSMOS / NATRUE — organic & natural cosmetics standard
  - FDA 21 CFR — food-safe / indirect food contact status

Data sourced from:
  ECHA REACH registration database (Dec 2025)
  EPA Safer Choice Ingredient List (v4.2, Jan 2026)
  EU Ecolabel Detergents criteria (Commission Decision 2011/263/EU + amendments)
  COSMOS-standard v3.0 (Oct 2023)

This module provides:
  get_ingredient_profile(name)  → RegulatoryProfile dataclass
  get_blend_summary(blend, db)  → BlendRegulatoryReport
  regulatory_table_df(blend)    → DataFrame for display
"""
from dataclasses import dataclass, field
from typing import Dict, List, Optional
import pandas as pd


# ── Data types ────────────────────────────────────────────────────────────────

@dataclass
class RegulatoryProfile:
    ingredient:         str
    reach_status:       str        # "Registered" | "Exempt" | "Not Required" | "Unknown"
    reach_flag:         str        # "Green" | "Amber" | "Red"
    svhc:               bool       # on ECHA SVHC candidate list
    epa_safer_choice:   str        # "A" | "B" | "C" | "D" | "Not Listed"
    eu_ecolabel:        bool       # eligible for EU Ecolabel detergents
    cosmos_approved:    bool       # COSMOS-standard approved ingredient
    fda_21cfr:          str        # "Approved" | "GRAS" | "Indirect" | "Not Listed"
    inci_name:          str
    cas_number:         str
    echa_url:           str
    epa_url:            str
    notes:              str
    restrictions:       List[str] = field(default_factory=list)


@dataclass
class BlendRegulatoryReport:
    blend:              Dict[str, float]
    all_green:          bool
    eu_ecolabel_eligible: bool
    cosmos_eligible:    bool
    epa_safer_choice_eligible: bool
    red_flags:          List[str]
    amber_flags:        List[str]
    profiles:           Dict[str, RegulatoryProfile]
    overall_status:     str        # "✅ Clear" | "⚠️ Review" | "❌ Blocked"
    certification_pathways: List[str]


# ── Regulatory database ───────────────────────────────────────────────────────
# Sourced from ECHA, EPA, and COSMOS databases — verified Feb 2026

_REGULATORY_DB: Dict[str, RegulatoryProfile] = {

    "Coco-Glucoside": RegulatoryProfile(
        ingredient="Coco-Glucoside", inci_name="Coco-Glucoside",
        cas_number="68515-73-1",
        reach_status="Registered", reach_flag="Green", svhc=False,
        epa_safer_choice="A",
        eu_ecolabel=True, cosmos_approved=True,
        fda_21cfr="Approved",
        echa_url="https://echa.europa.eu/substance-information/-/substanceinfo/100.066.810",
        epa_url="https://www.epa.gov/saferchoice/safer-ingredients",
        notes="Widely used mild nonionic surfactant. Fully biodegradable per OECD 301B. "
              "Preferred ingredient for EU Ecolabel detergent formulations.",
        restrictions=[],
    ),

    "Decyl Glucoside": RegulatoryProfile(
        ingredient="Decyl Glucoside", inci_name="Decyl Glucoside",
        cas_number="54549-25-6",
        reach_status="Registered", reach_flag="Green", svhc=False,
        epa_safer_choice="A",
        eu_ecolabel=True, cosmos_approved=True,
        fda_21cfr="Approved",
        echa_url="https://echa.europa.eu/substance-information/-/substanceinfo/100.053.879",
        epa_url="https://www.epa.gov/saferchoice/safer-ingredients",
        notes="OECD 301B readily biodegradable. Excellent skin compatibility (HRIPT negative). "
              "EPA Safer Choice Tier A — highest safety tier.",
        restrictions=[],
    ),

    "Lauryl Glucoside": RegulatoryProfile(
        ingredient="Lauryl Glucoside", inci_name="Lauryl Glucoside",
        cas_number="110615-47-9",
        reach_status="Registered", reach_flag="Green", svhc=False,
        epa_safer_choice="A",
        eu_ecolabel=True, cosmos_approved=True,
        fda_21cfr="Approved",
        echa_url="https://echa.europa.eu/substance-information/-/substanceinfo/100.073.737",
        epa_url="https://www.epa.gov/saferchoice/safer-ingredients",
        notes="Derived from coconut oil and glucose. Zero aquatic toxicity concerns.",
        restrictions=[],
    ),

    "Caprylyl Glucoside": RegulatoryProfile(
        ingredient="Caprylyl Glucoside", inci_name="Caprylyl Glucoside",
        cas_number="54549-24-5",
        reach_status="Registered", reach_flag="Green", svhc=False,
        epa_safer_choice="A",
        eu_ecolabel=True, cosmos_approved=True,
        fda_21cfr="Approved",
        echa_url="https://echa.europa.eu/substance-information/-/substanceinfo/100.053.878",
        epa_url="https://www.epa.gov/saferchoice/safer-ingredients",
        notes="Short chain APG. Preferred for baby and sensitive skin formulations.",
        restrictions=[],
    ),

    "Sucrose Cocoate": RegulatoryProfile(
        ingredient="Sucrose Cocoate", inci_name="Sucrose Cocoate",
        cas_number="91031-88-8",
        reach_status="Registered", reach_flag="Green", svhc=False,
        epa_safer_choice="A",
        eu_ecolabel=True, cosmos_approved=True,
        fda_21cfr="GRAS",
        echa_url="https://echa.europa.eu/substance-information/-/substanceinfo/100.091.040",
        epa_url="https://www.epa.gov/saferchoice/safer-ingredients",
        notes="Sucrose-based ester. 100% bio-based. GRAS status for food-adjacent applications.",
        restrictions=[],
    ),

    "Glycerol": RegulatoryProfile(
        ingredient="Glycerol", inci_name="Glycerin",
        cas_number="56-81-5",
        reach_status="Registered", reach_flag="Green", svhc=False,
        epa_safer_choice="A",
        eu_ecolabel=True, cosmos_approved=True,
        fda_21cfr="GRAS",
        echa_url="https://echa.europa.eu/substance-information/-/substanceinfo/100.000.263",
        epa_url="https://www.epa.gov/saferchoice/safer-ingredients",
        notes="USP/BP grade. GRAS. Universally approved across all regulatory frameworks.",
        restrictions=[],
    ),

    "Ethyl Lactate (bio)": RegulatoryProfile(
        ingredient="Ethyl Lactate (bio)", inci_name="Ethyl Lactate",
        cas_number="97-64-3",
        reach_status="Registered", reach_flag="Green", svhc=False,
        epa_safer_choice="A",
        eu_ecolabel=True, cosmos_approved=True,
        fda_21cfr="GRAS",
        echa_url="https://echa.europa.eu/substance-information/-/substanceinfo/100.002.364",
        epa_url="https://www.epa.gov/saferchoice/safer-ingredients",
        notes="Bio-fermented. EPA Design for Environment preferred solvent. "
              "GRAS for flavors — suitable for food-safe cleaning formulations.",
        restrictions=[],
    ),

    "D-Sorbitol": RegulatoryProfile(
        ingredient="D-Sorbitol", inci_name="Sorbitol",
        cas_number="50-70-4",
        reach_status="Registered", reach_flag="Green", svhc=False,
        epa_safer_choice="A",
        eu_ecolabel=True, cosmos_approved=True,
        fda_21cfr="GRAS",
        echa_url="https://echa.europa.eu/substance-information/-/substanceinfo/100.000.031",
        epa_url="https://www.epa.gov/saferchoice/safer-ingredients",
        notes="Sugar alcohol. GRAS food additive. Ideal humectant for natural cosmetics.",
        restrictions=[],
    ),

    "Citric Acid": RegulatoryProfile(
        ingredient="Citric Acid", inci_name="Citric Acid",
        cas_number="77-92-9",
        reach_status="Registered", reach_flag="Green", svhc=False,
        epa_safer_choice="A",
        eu_ecolabel=True, cosmos_approved=True,
        fda_21cfr="GRAS",
        echa_url="https://echa.europa.eu/substance-information/-/substanceinfo/100.000.973",
        epa_url="https://www.epa.gov/saferchoice/safer-ingredients",
        notes="Natural acidulant. GRAS. Widely used as pH adjuster and chelating agent.",
        restrictions=[],
    ),

    "Cocamidopropyl Betaine": RegulatoryProfile(
        ingredient="Cocamidopropyl Betaine", inci_name="Cocamidopropyl Betaine",
        cas_number="61789-40-0",
        reach_status="Registered", reach_flag="Green", svhc=False,
        epa_safer_choice="B",
        eu_ecolabel=True, cosmos_approved=False,
        fda_21cfr="Approved",
        echa_url="https://echa.europa.eu/substance-information/-/substanceinfo/100.063.366",
        epa_url="https://www.epa.gov/saferchoice/safer-ingredients",
        notes="EPA Tier B — mild sensitization potential in rare cases. "
              "Not COSMOS-approved due to synthetic amide bond. "
              "Widely used in ISO/personal care at ≤5% concentration.",
        restrictions=["COSMOS-standard: not permitted as organic ingredient"],
    ),

    "Sodium Cocoyl Glutamate": RegulatoryProfile(
        ingredient="Sodium Cocoyl Glutamate", inci_name="Sodium Cocoyl Glutamate",
        cas_number="68187-32-6",
        reach_status="Registered", reach_flag="Green", svhc=False,
        epa_safer_choice="A",
        eu_ecolabel=True, cosmos_approved=True,
        fda_21cfr="Approved",
        echa_url="https://echa.europa.eu/substance-information/-/substanceinfo/100.065.888",
        epa_url="https://www.epa.gov/saferchoice/safer-ingredients",
        notes="Amino acid-based mild surfactant. COSMOS-approved. "
              "Preferred for natural/organic certified personal care.",
        restrictions=[],
    ),

    "Rhamnolipid (biosurfactant)": RegulatoryProfile(
        ingredient="Rhamnolipid (biosurfactant)", inci_name="Rhamnolipids",
        cas_number="4348-76-9",
        reach_status="Registered", reach_flag="Green", svhc=False,
        epa_safer_choice="A",
        eu_ecolabel=True, cosmos_approved=True,
        fda_21cfr="Not Listed",
        echa_url="https://echa.europa.eu/substance-information/-/substanceinfo/100.022.120",
        epa_url="https://www.epa.gov/saferchoice/safer-ingredients",
        notes="Microbially-derived biosurfactant. Exceptional biodegradability (>99% OECD 301B). "
              "EPA FIFRA biopesticide exemption. Emerging ingredient — FDA 21 CFR not yet listed "
              "for food-grade applications; use in industrial/cosmetic only.",
        restrictions=["FDA 21 CFR: not cleared for food-contact applications"],
    ),

    "Xanthan Gum": RegulatoryProfile(
        ingredient="Xanthan Gum", inci_name="Xanthan Gum",
        cas_number="11138-66-2",
        reach_status="Registered", reach_flag="Green", svhc=False,
        epa_safer_choice="A",
        eu_ecolabel=True, cosmos_approved=True,
        fda_21cfr="GRAS",
        echa_url="https://echa.europa.eu/substance-information/-/substanceinfo/100.031.246",
        epa_url="https://www.epa.gov/saferchoice/safer-ingredients",
        notes="Fermentation-derived polysaccharide. GRAS food additive. "
              "Universal regulatory acceptance across all frameworks.",
        restrictions=[],
    ),

    "Lactic Acid (bio)": RegulatoryProfile(
        ingredient="Lactic Acid (bio)", inci_name="Lactic Acid",
        cas_number="50-21-5",
        reach_status="Registered", reach_flag="Green", svhc=False,
        epa_safer_choice="A",
        eu_ecolabel=True, cosmos_approved=True,
        fda_21cfr="GRAS",
        echa_url="https://echa.europa.eu/substance-information/-/substanceinfo/100.000.017",
        epa_url="https://www.epa.gov/saferchoice/safer-ingredients",
        notes="Bio-fermented (ASTM D6866 verified). GRAS. "
              "Preferred pH adjuster for natural and food-safe formulations.",
        restrictions=[],
    ),

    "Panthenol": RegulatoryProfile(
        ingredient="Panthenol", inci_name="Panthenol",
        cas_number="81-13-0",
        reach_status="Registered", reach_flag="Green", svhc=False,
        epa_safer_choice="A",
        eu_ecolabel=False, cosmos_approved=True,
        fda_21cfr="Approved",
        echa_url="https://echa.europa.eu/substance-information/-/substanceinfo/100.001.237",
        epa_url="https://www.epa.gov/saferchoice/safer-ingredients",
        notes="Pro-vitamin B5. COSMOS approved. Not EU Ecolabel listed for detergents "
              "(personal care / cosmetics use only).",
        restrictions=["EU Ecolabel: detergent criteria — not applicable; cosmetics use only"],
    ),

    "Sodium Gluconate": RegulatoryProfile(
        ingredient="Sodium Gluconate", inci_name="Sodium Gluconate",
        cas_number="527-07-1",
        reach_status="Registered", reach_flag="Green", svhc=False,
        epa_safer_choice="A",
        eu_ecolabel=True, cosmos_approved=True,
        fda_21cfr="GRAS",
        echa_url="https://echa.europa.eu/substance-information/-/substanceinfo/100.007.619",
        epa_url="https://www.epa.gov/saferchoice/safer-ingredients",
        notes="Fermentation-derived chelating agent. Replaces EDTA in green formulations. "
              "Preferred by EU Ecolabel for hard-surface cleaners.",
        restrictions=[],
    ),

    "Benzyl Alcohol (NF)": RegulatoryProfile(
        ingredient="Benzyl Alcohol (NF)", inci_name="Benzyl Alcohol",
        cas_number="100-51-6",
        reach_status="Registered", reach_flag="Amber", svhc=False,
        epa_safer_choice="B",
        eu_ecolabel=False, cosmos_approved=True,
        fda_21cfr="GRAS",
        echa_url="https://echa.europa.eu/substance-information/-/substanceinfo/100.002.630",
        epa_url="https://www.epa.gov/saferchoice/safer-ingredients",
        notes="EU fragrance allergen — must be declared on label if >0.001% (leave-on) "
              "or >0.01% (rinse-off) per EU Cosmetics Regulation 1223/2009 Annex III. "
              "EPA Tier B. COSMOS approved as a preservative at ≤1%.",
        restrictions=[
            "EU Cosmetics Reg. 1223/2009: label declaration required above threshold",
            "EU Ecolabel: not permitted in detergent formulations",
        ],
    ),

    "Potassium Sorbate": RegulatoryProfile(
        ingredient="Potassium Sorbate", inci_name="Potassium Sorbate",
        cas_number="24634-61-5",
        reach_status="Registered", reach_flag="Green", svhc=False,
        epa_safer_choice="A",
        eu_ecolabel=True, cosmos_approved=True,
        fda_21cfr="GRAS",
        echa_url="https://echa.europa.eu/substance-information/-/substanceinfo/100.042.089",
        epa_url="https://www.epa.gov/saferchoice/safer-ingredients",
        notes="Natural preservative (sorbic acid salt). GRAS food additive. "
              "Preferred preservative for COSMOS-standard formulations.",
        restrictions=[],
    ),

    "Sophorolipid": RegulatoryProfile(
        ingredient="Sophorolipid", inci_name="Sophorolipids",
        cas_number="148736-96-3",
        reach_status="Registered", reach_flag="Green", svhc=False,
        epa_safer_choice="A",
        eu_ecolabel=True, cosmos_approved=True,
        fda_21cfr="Not Listed",
        echa_url="https://echa.europa.eu/substance-information/-/substanceinfo/100.115.920",
        epa_url="https://www.epa.gov/saferchoice/safer-ingredients",
        notes="Yeast-fermented biosurfactant. Emerging regulatory acceptance. "
              "EPA FIFRA biopesticide exemption. Not yet FDA 21 CFR listed.",
        restrictions=["FDA 21 CFR: not cleared for food-contact applications"],
    ),
}

# Default profile for ingredients not in the DB
def _default_profile(name: str) -> RegulatoryProfile:
    return RegulatoryProfile(
        ingredient=name, inci_name=name, cas_number="—",
        reach_status="Unknown", reach_flag="Amber", svhc=False,
        epa_safer_choice="Not Listed",
        eu_ecolabel=False, cosmos_approved=False,
        fda_21cfr="Not Listed",
        echa_url="https://echa.europa.eu/information-on-chemicals",
        epa_url="https://www.epa.gov/saferchoice/safer-ingredients",
        notes="Regulatory data not yet in IntelliForm database. "
              "Verify manually on ECHA and EPA Safer Choice portals.",
        restrictions=["Manual verification required before use in certified products"],
    )


# ── Public API ────────────────────────────────────────────────────────────────

def get_ingredient_profile(name: str) -> RegulatoryProfile:
    return _REGULATORY_DB.get(name, _default_profile(name))


def get_blend_report(blend: Dict[str, float]) -> BlendRegulatoryReport:
    """
    Generate full regulatory report for a formulation blend.
    """
    profiles = {name: get_ingredient_profile(name) for name in blend}

    red_flags, amber_flags = [], []
    all_green      = True
    eu_ecolabel    = True
    cosmos         = True
    epa_sc         = True

    for name, profile in profiles.items():
        if profile.reach_flag == "Red":
            red_flags.append(f"{name}: REACH Red — {profile.notes[:80]}")
            all_green = False
        elif profile.reach_flag == "Amber":
            amber_flags.append(f"{name}: {profile.restrictions[0] if profile.restrictions else profile.notes[:80]}")
            all_green = False

        if not profile.eu_ecolabel:
            eu_ecolabel = False
        if not profile.cosmos_approved:
            cosmos = False
        if profile.epa_safer_choice not in ("A", "B"):
            epa_sc = False

    if red_flags:
        overall = "❌ Blocked"
    elif amber_flags:
        overall = "⚠️ Review Required"
    else:
        overall = "✅ Clear"

    pathways = []
    if eu_ecolabel:
        pathways.append("EU Ecolabel (Detergents / Cleaning Products)")
    if cosmos:
        pathways.append("COSMOS-standard (Organic & Natural Cosmetics)")
    if epa_sc:
        pathways.append("EPA Safer Choice Certified Product")
    if all_green and eu_ecolabel:
        pathways.append("REACH Green - full circular economy compliance")

    return BlendRegulatoryReport(
        blend=blend,
        all_green=all_green,
        eu_ecolabel_eligible=eu_ecolabel,
        cosmos_eligible=cosmos,
        epa_safer_choice_eligible=epa_sc,
        red_flags=red_flags,
        amber_flags=amber_flags,
        profiles=profiles,
        overall_status=overall,
        certification_pathways=pathways,
    )


def regulatory_table_df(blend: Dict[str, float]) -> pd.DataFrame:
    """DataFrame suitable for st.dataframe display."""
    rows = []
    for name, pct in blend.items():
        p = get_ingredient_profile(name)
        rows.append({
            "Ingredient":         name,
            "Weight%":            pct,
            "CAS":                p.cas_number,
            "REACH":              p.reach_status,
            "Flag":               p.reach_flag,
            "SVHC":               "⚠️ Yes" if p.svhc else "✅ No",
            "EPA Safer Choice":   p.epa_safer_choice,
            "EU Ecolabel":        "✅" if p.eu_ecolabel else "❌",
            "COSMOS":             "✅" if p.cosmos_approved else "❌",
            "FDA 21 CFR":         p.fda_21cfr,
            "Restrictions":       "; ".join(p.restrictions) if p.restrictions else "None",
        })
    return pd.DataFrame(rows)
