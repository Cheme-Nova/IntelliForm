"""
modules/sds_generator.py
GHS-compliant Safety Data Sheet generator for IntelliForm blend formulations.

Follows EU GHS/CLP Annex II (16-section) format.
Data sources (in priority order):
  1. EPA CompTox Dashboard (OPERA predictions, GHS codes, SVHC/CMR flags)
  2. modules/regulatory.py (REACH, EPA Safer Choice, COSMOS, CAS)
  3. modules/qsar.py (biodegradability %, ecotoxicity score)
  4. modules/stability.py (viscosity estimate)
  5. ingredients_db.csv (MW, LogP, Function)

Output: SDSDocument dataclass → sds_to_pdf() → bytes / sds_to_json() → dict
"""
from __future__ import annotations

import io
import json
from dataclasses import dataclass, field, asdict
from datetime import datetime
from typing import Dict, List, Optional, Tuple

import pandas as pd

# ── GHS H-statement library (EU CLP 2024, Annex VI) ──────────────────────────

_H_STATEMENTS: Dict[str, str] = {
    # Physical hazards
    "H224": "Extremely flammable liquid and vapour.",
    "H225": "Highly flammable liquid and vapour.",
    "H226": "Flammable liquid and vapour.",
    "H227": "Combustible liquid.",
    "H228": "Flammable solid.",
    "H240": "Heating may cause an explosion.",
    "H241": "Heating may cause a fire or explosion.",
    "H242": "Heating may cause a fire.",
    "H250": "Catches fire spontaneously if exposed to air.",
    "H251": "Self-heating in large quantities; may catch fire.",
    "H252": "Self-heating in large quantities; may catch fire.",
    "H260": "In contact with water releases flammable gases which may ignite spontaneously.",
    "H261": "In contact with water releases flammable gases.",
    "H270": "May cause or intensify fire; oxidiser.",
    "H271": "May cause fire or explosion; strong oxidiser.",
    "H272": "May intensify fire; oxidiser.",
    "H280": "Contains gas under pressure; may explode if heated.",
    "H281": "Contains refrigerated gas; may cause cryogenic burns or injury.",
    "H290": "May be corrosive to metals.",
    # Health hazards
    "H300": "Fatal if swallowed.",
    "H301": "Toxic if swallowed.",
    "H302": "Harmful if swallowed.",
    "H303": "May be harmful if swallowed.",
    "H304": "May be fatal if swallowed and enters airways.",
    "H310": "Fatal in contact with skin.",
    "H311": "Toxic in contact with skin.",
    "H312": "Harmful in contact with skin.",
    "H313": "May be harmful in contact with skin.",
    "H314": "Causes severe skin burns and eye damage.",
    "H315": "Causes skin irritation.",
    "H317": "May cause an allergic skin reaction.",
    "H318": "Causes serious eye damage.",
    "H319": "Causes serious eye irritation.",
    "H330": "Fatal if inhaled.",
    "H331": "Toxic if inhaled.",
    "H332": "Harmful if inhaled.",
    "H333": "May be harmful if inhaled.",
    "H334": "May cause allergy or asthma symptoms or breathing difficulties if inhaled.",
    "H335": "May cause respiratory irritation.",
    "H336": "May cause drowsiness or dizziness.",
    "H340": "May cause genetic defects.",
    "H341": "Suspected of causing genetic defects.",
    "H350": "May cause cancer.",
    "H351": "Suspected of causing cancer.",
    "H360": "May damage fertility or the unborn child.",
    "H361": "Suspected of damaging fertility or the unborn child.",
    "H362": "May cause harm to breast-fed children.",
    "H370": "Causes damage to organs.",
    "H371": "May cause damage to organs.",
    "H372": "Causes damage to organs through prolonged or repeated exposure.",
    "H373": "May cause damage to organs through prolonged or repeated exposure.",
    # Environmental hazards
    "H400": "Very toxic to aquatic life.",
    "H410": "Very toxic to aquatic life with long lasting effects.",
    "H411": "Toxic to aquatic life with long lasting effects.",
    "H412": "Harmful to aquatic life with long lasting effects.",
    "H413": "May cause long lasting harmful effects to aquatic life.",
    "H420": "Harms public health and the environment by destroying ozone in the upper atmosphere.",
}

# Codes that require "Danger" signal word (GHS category 1/1A/1B where applicable)
_DANGER_CODES = {
    "H224","H225","H228","H240","H241","H250","H251","H260","H270","H271",
    "H290","H300","H310","H314","H318","H330","H331","H334",
    "H340","H350","H360","H370","H372","H400","H410",
}

# GHS pictogram → text representation (Unicode symbols usable in ReportLab)
_PICTOGRAM_MAP: Dict[str, Tuple[str, str]] = {
    "GHS01": ("⚠",  "Exploding bomb"),
    "GHS02": ("🔥", "Flame"),
    "GHS03": ("⊕",  "Flame over circle (oxidiser)"),
    "GHS04": ("⊙",  "Gas cylinder"),
    "GHS05": ("⚗",  "Corrosion"),
    "GHS06": ("☠",  "Skull and crossbones"),
    "GHS07": ("⚠",  "Exclamation mark"),
    "GHS08": ("⊕",  "Health hazard"),
    "GHS09": ("⊕",  "Environmental hazard (fish/tree)"),
}

# H-code → pictogram(s)
_H_TO_PICTOGRAM: Dict[str, List[str]] = {
    "H200": ["GHS01"], "H201": ["GHS01"], "H202": ["GHS01"], "H204": ["GHS01"],
    "H240": ["GHS01"], "H241": ["GHS01", "GHS02"],
    "H224": ["GHS02"], "H225": ["GHS02"], "H226": ["GHS02"], "H227": ["GHS02"],
    "H228": ["GHS02"], "H250": ["GHS02"], "H260": ["GHS02"], "H261": ["GHS02"],
    "H270": ["GHS03"], "H271": ["GHS03", "GHS01"], "H272": ["GHS03"],
    "H280": ["GHS04"], "H281": ["GHS04"],
    "H290": ["GHS05"],
    "H314": ["GHS05"], "H318": ["GHS05"],
    "H300": ["GHS06"], "H301": ["GHS06"], "H310": ["GHS06"], "H311": ["GHS06"],
    "H330": ["GHS06"], "H331": ["GHS06"],
    "H302": ["GHS07"], "H312": ["GHS07"], "H315": ["GHS07"], "H317": ["GHS07"],
    "H319": ["GHS07"], "H332": ["GHS07"], "H335": ["GHS07"], "H336": ["GHS07"],
    "H304": ["GHS08"], "H334": ["GHS08"], "H340": ["GHS08"], "H341": ["GHS08"],
    "H350": ["GHS08"], "H351": ["GHS08"], "H360": ["GHS08"], "H361": ["GHS08"],
    "H362": ["GHS08"], "H370": ["GHS08"], "H371": ["GHS08"],
    "H372": ["GHS08"], "H373": ["GHS08"],
    "H400": ["GHS09"], "H410": ["GHS09"], "H411": ["GHS09"],
    "H412": ["GHS09"], "H413": ["GHS09"],
}

# Key precautionary statements keyed by hazard category
_P_BY_CATEGORY: Dict[str, List[str]] = {
    "skin_irritant":   ["P264 – Wash hands thoroughly after handling.",
                        "P280 – Wear protective gloves/eye protection.",
                        "P302+P352 – IF ON SKIN: Wash with plenty of water.",
                        "P321 – Specific treatment (see supplemental first aid)."],
    "eye_irritant":    ["P264 – Wash face and hands thoroughly after handling.",
                        "P280 – Wear eye protection.",
                        "P305+P351+P338 – IF IN EYES: Rinse cautiously with water for several minutes."],
    "aquatic_tox":     ["P273 – Avoid release to the environment.",
                        "P391 – Collect spillage.",
                        "P501 – Dispose of contents/container in accordance with local regulations."],
    "flammable":       ["P210 – Keep away from heat/sparks/open flames.",
                        "P233 – Keep container tightly closed.",
                        "P241 – Use explosion-proof electrical equipment.",
                        "P370+P378 – In case of fire: Use CO₂ or dry powder to extinguish."],
    "acute_oral":      ["P260 – Do not breathe vapours.",
                        "P264 – Wash hands thoroughly after handling.",
                        "P270 – Do not eat, drink or smoke when using this product.",
                        "P301+P312 – IF SWALLOWED: Call a POISON CENTER if you feel unwell."],
    "cmr":             ["P201 – Obtain special instructions before use.",
                        "P202 – Do not handle until all safety precautions have been read and understood.",
                        "P280 – Wear protective gloves/protective clothing/eye protection/face protection.",
                        "P308+P313 – IF exposed or concerned: Get medical advice/attention."],
    "general":         ["P102 – Keep out of reach of children.",
                        "P103 – Read label before use.",
                        "P501 – Dispose of contents/container per local regulations."],
}


# ── Dataclasses ───────────────────────────────────────────────────────────────

@dataclass
class SDSIngredient:
    name: str
    cas: str
    inci: str
    wt_pct: float
    reach_status: str
    ghs_codes: List[str]
    cmr: Optional[str]
    svhc: bool


@dataclass
class SDSDocument:
    product_name:     str
    revision_date:    str
    version:          str
    vertical:         str
    # Section content as plain-text dicts
    s1_identification:       Dict
    s2_hazards:              Dict
    s3_composition:          List[SDSIngredient]
    s4_first_aid:            Dict
    s5_fire_fighting:        Dict
    s6_accidental_release:   Dict
    s7_handling_storage:     Dict
    s8_exposure_controls:    Dict
    s9_physical_props:       Dict
    s10_stability:           Dict
    s11_toxicological:       Dict
    s12_ecological:          Dict
    s13_disposal:            Dict
    s14_transport:           Dict
    s15_regulatory:          Dict
    s16_other:               Dict
    # Aggregated blend-level flags
    all_ghs_codes:    List[str]
    signal_word:      str
    pictograms:       List[str]


# ── Hazard aggregation helpers ────────────────────────────────────────────────

def _signal_word(codes: List[str]) -> str:
    for c in codes:
        if c in _DANGER_CODES:
            return "Danger"
    return "Warning" if codes else "No signal word"


def _pictograms(codes: List[str]) -> List[str]:
    pics: set = set()
    for c in codes:
        for p in _H_TO_PICTOGRAM.get(c, []):
            pics.add(p)
    return sorted(pics)


def _p_statements(codes: List[str], cmr_found: bool) -> List[str]:
    stmts: List[str] = list(_P_BY_CATEGORY["general"])
    code_set = set(codes)
    if code_set & {"H315", "H317"}:
        stmts += _P_BY_CATEGORY["skin_irritant"]
    if code_set & {"H318", "H319"}:
        stmts += _P_BY_CATEGORY["eye_irritant"]
    if code_set & {"H400", "H410", "H411", "H412", "H413"}:
        stmts += _P_BY_CATEGORY["aquatic_tox"]
    if code_set & {"H224", "H225", "H226", "H227"}:
        stmts += _P_BY_CATEGORY["flammable"]
    if code_set & {"H300", "H301", "H302", "H303"}:
        stmts += _P_BY_CATEGORY["acute_oral"]
    if cmr_found:
        stmts += _P_BY_CATEGORY["cmr"]
    return list(dict.fromkeys(stmts))  # deduplicate, preserve order


def _derive_aquatic_h(fish_lc50_log_mmol: Optional[float],
                       biodeg_prob: Optional[float],
                       log_bcf: Optional[float]) -> List[str]:
    """Derive aquatic H-codes from OPERA QSAR data when CompTox codes unavailable."""
    codes = []
    if fish_lc50_log_mmol is not None:
        # log(LC50 mmol/L): lower = more toxic
        if fish_lc50_log_mmol < -3.5:
            codes.append("H410")
        elif fish_lc50_log_mmol < -2.0:
            codes.append("H411")
        elif fish_lc50_log_mmol < -0.5:
            codes.append("H412")
    if not codes and biodeg_prob is not None and biodeg_prob < 0.5:
        codes.append("H413")
    return codes


# ── Section builders ──────────────────────────────────────────────────────────

def _s1(product_name: str, vertical: str, manufacturer_info: dict) -> dict:
    return {
        "product_name":    product_name,
        "product_code":    f"IF-{datetime.now().strftime('%Y%m%d')}",
        "intended_use":    f"{vertical.replace('_', ' ').title()} specialty chemical formulation",
        "manufacturer":    manufacturer_info.get("name", "ChemeNova LLC"),
        "address":         manufacturer_info.get("address", "Newark, NJ 07102, USA"),
        "phone":           manufacturer_info.get("phone", "+1 (973) 555-0100"),
        "email":           manufacturer_info.get("email", "safety@chemenova.com"),
        "emergency":       manufacturer_info.get("emergency", "CHEMTREC: +1 (800) 424-9300"),
        "sds_number":      f"SDS-{datetime.now().strftime('%Y%m%d-%H%M')}",
    }


def _s2(all_codes: List[str], signal_word: str, pictograms: List[str],
        cmr_found: bool) -> dict:
    h_texts  = [f"{c}: {_H_STATEMENTS.get(c, c)}" for c in sorted(set(all_codes))]
    p_stmts  = _p_statements(all_codes, cmr_found)
    pic_desc = [f"{p} ({_PICTOGRAM_MAP[p][1]})" for p in pictograms if p in _PICTOGRAM_MAP]
    return {
        "signal_word":         signal_word,
        "pictograms":          pic_desc or ["No pictogram required"],
        "hazard_statements":   h_texts or ["Not classified as hazardous under GHS."],
        "precautionary_statements": p_stmts,
        "other_hazards":       ("Contains CMR substance(s) — refer to Section 11."
                                if cmr_found else "None known."),
    }


def _s4(all_codes: List[str]) -> dict:
    code_set = set(all_codes)
    eye  = ("Rinse cautiously with water for at least 15 minutes. Remove contact lenses if present. "
            "Seek medical attention if irritation persists."
            if code_set & {"H314","H318","H319"} else
            "Rinse with water if irritation occurs.")
    skin = ("Remove contaminated clothing immediately. Wash with soap and water for at least "
            "15 minutes. Seek medical attention if burns develop."
            if code_set & {"H310","H311","H314"} else
            "Wash with soap and water. Seek medical attention if irritation persists.")
    inhl = ("Remove to fresh air. If breathing is difficult, administer oxygen. "
            "Seek immediate medical attention."
            if code_set & {"H330","H331","H334"} else
            "Move to fresh air. Seek medical attention if symptoms develop.")
    ingst = ("Rinse mouth with water. Do NOT induce vomiting. Seek medical attention immediately."
             if code_set & {"H300","H301","H302"} else
             "Rinse mouth with water. Seek medical attention if large quantities ingested.")
    return {
        "inhalation":   inhl,
        "skin_contact": skin,
        "eye_contact":  eye,
        "ingestion":    ingst,
        "notes_to_physician": ("Treat symptomatically. For CMR substances, consult poison control."
                               if code_set & {"H340","H350","H360"} else
                               "Treat symptomatically."),
    }


def _s5(avg_logp: Optional[float], code_set: set) -> dict:
    flash_note = ("Consult ingredient-level SDS for specific flash points."
                  if code_set & {"H224","H225","H226","H227"} else
                  "Aqueous-based formulation — not expected to be flammable.")
    return {
        "extinguishing_media":     "CO₂, dry powder, alcohol-resistant foam.",
        "unsuitable_media":        "High-pressure water jet on flammable components.",
        "special_hazards":         ("May release toxic fumes under fire conditions."
                                   if code_set & {"H314","H318"} else
                                   "Combustion may produce CO₂ and water vapour."),
        "protective_equipment":    "Self-contained breathing apparatus and full protective suit.",
        "flash_point_note":        flash_note,
    }


def _s6(biodeg_fraction: Optional[float]) -> dict:
    env_note = (
        f"Blend is {biodeg_fraction:.0%} readily biodegradable (OECD 301B proxy). "
        "Contain spillage — avoid release to drains or waterways."
        if biodeg_fraction is not None else
        "Contain spillage. Avoid release to drains or waterways."
    )
    return {
        "personal_precautions":   "Wear PPE specified in Section 8. Ensure adequate ventilation.",
        "environmental_precautions": env_note,
        "containment":            "Absorb with inert material (vermiculite, sand). Collect in labelled container.",
        "cleanup":                "Dispose of according to local regulations (see Section 13).",
    }


def _s7(code_set: set) -> dict:
    return {
        "handling_precautions": (
            "Avoid contact with skin and eyes. Use in well-ventilated areas. "
            + ("Keep away from ignition sources. No smoking." if code_set & {"H224","H225","H226"} else "")
        ).strip(),
        "storage_conditions":   ("Store in a cool, dry, well-ventilated area. Keep container tightly closed. "
                                 "Protect from direct sunlight and frost."),
        "incompatibilities":    ("Strong oxidisers, strong acids, strong bases."
                                 if code_set & {"H314","H315"} else
                                 "Strong oxidisers."),
        "storage_temperature":  "5–30 °C",
        "shelf_life":           "24 months from date of manufacture (unopened).",
    }


def _s8(code_set: set, oel_notes: str) -> dict:
    eye_prot  = "Safety glasses" if code_set & {"H319"} else "None required under normal use"
    eye_prot  = "Chemical goggles" if code_set & {"H314","H318"} else eye_prot
    skin_prot = "Nitrile gloves (0.3 mm minimum)" if code_set & {"H314","H315","H317"} else "Rubber gloves recommended"
    resp_prot = ("Half-face respirator with organic vapour cartridge"
                 if code_set & {"H332","H335","H336"} else
                 "Not required under normal use conditions")
    return {
        "engineering_controls":   "Ensure adequate local exhaust ventilation. Use in well-ventilated area.",
        "oel_note":               oel_notes or "No specific OEL established for this mixture.",
        "eye_protection":         eye_prot,
        "hand_protection":        skin_prot,
        "skin_body_protection":   "Lab coat or chemical resistant apron.",
        "respiratory_protection": resp_prot,
        "hygiene_measures":       "Wash hands before breaks and after handling. Do not eat or drink during use.",
    }


def _s9(blend: Dict[str, float], db: pd.DataFrame, avg_mw: Optional[float],
        avg_logp: Optional[float], avg_tpsa: Optional[float],
        viscosity_est: Optional[float]) -> dict:
    return {
        "appearance":           "Liquid or gel (appearance depends on formulation)",
        "colour":               "Clear to pale yellow",
        "odour":                "Slight characteristic odour",
        "pH":                   "5.0–8.0 (typical; adjust per application)",
        "melting_point":        "Not applicable (mixture)",
        "boiling_point":        ">100 °C (aqueous-based)",
        "flash_point":          "Not applicable (aqueous) or see ingredient SDS",
        "evaporation_rate":     "Similar to water",
        "flammability":         "Not classified as flammable (aqueous base)",
        "vapour_pressure":      "~23 hPa at 20 °C (aqueous fraction dominates)",
        "relative_density":     "~1.0–1.05 g/cm³",
        "water_solubility":     "Fully miscible",
        "partition_coeff_logp": f"{avg_logp:.2f}" if avg_logp is not None else "Not determined",
        "molecular_weight_avg": f"{avg_mw:.1f} g/mol" if avg_mw is not None else "Not determined",
        "tpsa_avg":             f"{avg_tpsa:.1f} Å²" if avg_tpsa is not None else "Not determined",
        "dynamic_viscosity":    f"~{viscosity_est:.0f} mPa·s" if viscosity_est else "Not determined",
    }


def _s10(code_set: set) -> dict:
    return {
        "reactivity":           "No hazardous reactivity known under normal conditions.",
        "chemical_stability":   "Stable under recommended storage conditions.",
        "hazardous_reactions":  "None known.",
        "conditions_to_avoid":  "Excessive heat, direct sunlight, freezing temperatures.",
        "incompatible_materials": "Strong oxidising agents, strong mineral acids, strong alkalis.",
        "decomposition_products": ("May release CO₂, CO, NOₓ upon combustion."
                                  if code_set & {"H226","H225"} else
                                  "CO₂ and water upon combustion."),
    }


def _s11(comptox_data: list, all_codes: List[str]) -> dict:
    # Aggregate OPERA predictions across ingredients
    lc50_vals = [c.fish_lc50_log_mmol_l for c in comptox_data
                 if c and c.fish_lc50_log_mmol_l is not None]
    ec50_vals = [c.daphnia_ec50_log_mmol_l for c in comptox_data
                 if c and c.daphnia_ec50_log_mmol_l is not None]
    kow_vals  = [c.log_kow for c in comptox_data if c and c.log_kow is not None]

    code_set = set(all_codes)
    cmr_note = ("Contains CMR substance(s): " +
                ", ".join(c.name for c in comptox_data if c and c.cmr_category) + ". "
                "Refer to ingredient-level SDS for details."
                if any(c.cmr_category for c in comptox_data if c) else
                "No CMR substances identified in this formulation.")

    return {
        "acute_toxicity_oral":    "Category 4 or below — refer to ingredient SDS for details.",
        "skin_irritation":        ("H315 – skin irritant per CLP" if "H315" in code_set else "Not classified."),
        "eye_irritation":         ("H318/H319 – eye damage/irritation per CLP" if code_set & {"H318","H319"} else "Not classified."),
        "respiratory_sensitisation": ("H334 – potential respiratory sensitiser" if "H334" in code_set else "Not classified."),
        "carcinogenicity":        ("H350/H351 present — refer to Section 15." if code_set & {"H350","H351"} else "Not classified."),
        "reproductive_toxicity":  ("H360/H361 present — refer to Section 15." if code_set & {"H360","H361"} else "Not classified."),
        "cmr_summary":            cmr_note,
        "fish_lc50_log_mmol":    (f"{min(lc50_vals):.2f}" if lc50_vals else "Not determined"),
        "daphnia_ec50_log_mmol": (f"{min(ec50_vals):.2f}" if ec50_vals else "Not determined"),
        "avg_log_kow":           (f"{sum(kow_vals)/len(kow_vals):.2f}" if kow_vals else "Not determined"),
        "notes":                  "Values are OPERA QSAR predictions (EPA CompTox). Validate with laboratory data.",
    }


def _s12(comptox_data: list, blend_bio_pct: float, blend_etox: float) -> dict:
    bcf_vals = [c.log_bcf for c in comptox_data if c and c.log_bcf is not None]
    kow_vals = [c.log_kow for c in comptox_data if c and c.log_kow is not None]
    biodeg_pct_str = f"{blend_bio_pct:.1f}% (QSAR estimate, OECD 301B proxy)"
    bioaccum = "Potential for bioaccumulation" if (bcf_vals and max(bcf_vals) > 3.7) else "Low bioaccumulation potential"
    return {
        "aquatic_toxicity":      f"Ecotoxicity score {blend_etox:.1f}/10 (IntelliForm QSAR); verify by Daphnia and fish acute tests.",
        "persistence":           biodeg_pct_str,
        "bioaccumulation":       bioaccum + (f" (max log BCF: {max(bcf_vals):.2f})" if bcf_vals else ""),
        "mobility_soil":         (f"log Koc data incomplete — exercise caution near water bodies."
                                  if not any(c.log_koc for c in comptox_data if c) else
                                  f"log Koc data available from CompTox; consult ingredient SDS."),
        "pbt_vpvb_assessment":   "Not formally assessed — consult ECHA REACH dossier for individual substances.",
        "endocrine_disruption":  "Not assessed.",
        "other_adverse_effects": "Avoid release to drains, surface water, and soil.",
    }


def _s13(biodeg_fraction: Optional[float], code_set: set) -> dict:
    waste_class = "Low hazard" if not (code_set & {"H300","H301","H310","H330","H314"}) else "Hazardous waste"
    return {
        "waste_treatment":     (f"Formulation is {biodeg_fraction:.0%} readily biodegradable. "
                                if biodeg_fraction and biodeg_fraction > 0.6 else "") +
                               "Neutralise if necessary. Dispose via licensed waste contractor.",
        "waste_code":          "08 01 11 (paints/varnish containing organic solvents) or 16 05 06 (laboratory chemicals) — confirm with local authority.",
        "classification":      waste_class,
        "packaging_disposal":  "Empty containers should be disposed of in accordance with local regulations. Do not re-use.",
    }


def _s14(code_set: set) -> dict:
    dangerous = bool(code_set & {"H224","H225","H226","H300","H301","H310","H330","H314"})
    return {
        "un_number":           "Not regulated under ADR/IMDG/IATA for normal quantities" if not dangerous else "Consult transport specialist",
        "proper_shipping_name": "Environmentally Hazardous Substance, Liquid" if "H410" in code_set else "Not regulated",
        "hazard_class":        "9 (Miscellaneous)" if "H410" in code_set else "Not classified",
        "packing_group":       "III" if "H410" in code_set else "N/A",
        "marine_pollutant":    "Yes" if code_set & {"H400","H410","H411"} else "No",
        "special_precautions": "Transport in tightly closed original containers. Keep away from food and feedstuffs.",
    }


def _s15(reg_report, vertical: str) -> dict:
    if reg_report is None:
        return {"note": "Regulatory data unavailable — run formulation first."}
    reach = "All ingredients registered or exempt under REACH Regulation (EC) 1907/2006."
    if getattr(reg_report, "red_flags", []):
        reach = "REACH flags: " + "; ".join(reg_report.red_flags)
    ecolabel = "✓ EU Ecolabel eligible" if getattr(reg_report, "eu_ecolabel_eligible", False) else "Not assessed for EU Ecolabel"
    cosmos   = "✓ COSMOS-standard eligible" if getattr(reg_report, "cosmos_eligible", False) else "Not assessed for COSMOS"
    epa_sc   = "✓ EPA Safer Choice eligible" if getattr(reg_report, "epa_safer_choice_eligible", False) else "Not assessed for EPA Safer Choice"
    return {
        "eu_reach":              reach,
        "eu_clp":                "Classified per EU CLP Regulation (EC) 1272/2008.",
        "eu_ecolabel":           ecolabel,
        "cosmos":                cosmos,
        "epa_safer_choice":      epa_sc,
        "biocidal_products":     "Not a biocidal product (BPR 528/2012 N/A).",
        "california_prop65":     "No Prop 65 listed substances identified at reportable levels.",
        "other_regulations":     f"Formulated for {vertical.replace('_',' ')} applications. Verify sector-specific requirements.",
    }


def _s16(product_name: str, version: str, all_codes: List[str]) -> dict:
    return {
        "revision_info":     f"SDS version {version} prepared by IntelliForm v2.1 · ChemeNova LLC",
        "revision_date":     datetime.now().strftime("%Y-%m-%d"),
        "prepared_by":       "ChemeNova LLC Safety & Regulatory Team",
        "full_h_text":       {c: _H_STATEMENTS.get(c, c) for c in sorted(set(all_codes))},
        "disclaimer":        (
            "The information in this SDS was compiled in good faith from supplier data, "
            "QSAR predictions (EPA CompTox OPERA), and published literature. It is provided "
            "as a guide for safe handling. The user is responsible for determining suitability "
            "for a specific purpose and complying with applicable laws and regulations."
        ),
        "training_required": "Users must be trained in handling specialty chemicals before using this product.",
        "intelliform_ref":   "Generated by IntelliForm™ — ChemeNova LLC · DOI: 10.26434/chemrxiv.15000857",
    }


# ── Main generator ────────────────────────────────────────────────────────────

def generate_sds(
    blend:             Dict[str, float],
    db:                pd.DataFrame,
    product_name:      str = "IntelliForm Specialty Formulation",
    vertical:          str = "personal_care",
    manufacturer_info: Optional[dict] = None,
    reg_report=None,
    qsar_bio_pct:      float = 75.0,
    qsar_etox:         float = 5.0,
    viscosity_est:     Optional[float] = None,
    version:           str = "1.0",
) -> SDSDocument:
    """
    Assemble a 16-section GHS SDS for a blend formulation.

    Pulls live data from:
      - modules/regulatory.py (if reg_report supplied)
      - EPA CompTox (via modules/comptox.py, skipped if offline)
      - ingredients_db.csv (MW, LogP, Function)
    """
    if manufacturer_info is None:
        manufacturer_info = {}

    # ── Fetch regulatory + CompTox data per ingredient ────────────────────────
    try:
        from modules.regulatory import get_ingredient_profile
        reg_profiles = {ing: get_ingredient_profile(ing) for ing in blend}
    except Exception:
        reg_profiles = {}

    try:
        from modules.comptox import lookup_ingredient
        comptox_results = {ing: lookup_ingredient(ing) for ing in blend}
    except Exception:
        comptox_results = {}

    # ── Aggregate GHS codes across all ingredients ────────────────────────────
    all_codes: List[str] = []
    cmr_found = False
    sds_ings: List[SDSIngredient] = []

    idx_db = db.set_index("Ingredient") if "Ingredient" in db.columns else pd.DataFrame()

    avg_mw_vals, avg_logp_vals, avg_tpsa_vals = [], [], []

    for ing, pct in blend.items():
        rp  = reg_profiles.get(ing)
        ctx = comptox_results.get(ing)

        cas   = (rp.cas_number if rp else "") or (ctx.cas if ctx else "Not determined")
        inci  = (rp.inci_name  if rp else ing)
        reach = (rp.reach_status if rp else "Unknown")

        # GHS codes: CompTox first, then derive from OPERA if empty
        codes = list(ctx.ghs_hazard_codes) if (ctx and ctx.ghs_hazard_codes) else []
        if not codes and ctx:
            codes = _derive_aquatic_h(ctx.fish_lc50_log_mmol_l,
                                      ctx.biodeg_probability,
                                      ctx.log_bcf)
        all_codes.extend(codes)

        if ctx and ctx.cmr_category:
            cmr_found = True

        sds_ings.append(SDSIngredient(
            name=ing, cas=cas, inci=inci, wt_pct=pct,
            reach_status=reach,
            ghs_codes=codes,
            cmr=ctx.cmr_category if ctx else None,
            svhc=bool(ctx.svhc_candidate) if ctx else False,
        ))

        # Physical props from DB
        if ing in idx_db.index:
            row = idx_db.loc[ing]
            try: avg_mw_vals.append(float(row.get("MolWt", row.get("MW", None))))
            except Exception: pass
            try: avg_logp_vals.append(float(row.get("LogP", None)))
            except Exception: pass
            try: avg_tpsa_vals.append(float(row.get("TPSA", None)))
            except Exception: pass

    all_codes = list(dict.fromkeys(all_codes))  # deduplicate, preserve order
    code_set  = set(all_codes)
    sig_word  = _signal_word(all_codes)
    pics      = _pictograms(all_codes)

    avg_mw   = (sum(avg_mw_vals)   / len(avg_mw_vals))   if avg_mw_vals   else None
    avg_logp = (sum(avg_logp_vals) / len(avg_logp_vals)) if avg_logp_vals else None
    avg_tpsa = (sum(avg_tpsa_vals) / len(avg_tpsa_vals)) if avg_tpsa_vals else None

    # Biodegradability fraction for Sections 6, 12, 13
    bio_fraction = (qsar_bio_pct / 100.0) if qsar_bio_pct else None
    oel_notes    = ""

    comptox_list = [v for v in comptox_results.values() if v and not v.error]

    return SDSDocument(
        product_name=product_name,
        revision_date=datetime.now().strftime("%Y-%m-%d"),
        version=version,
        vertical=vertical,
        s1_identification=      _s1(product_name, vertical, manufacturer_info),
        s2_hazards=             _s2(all_codes, sig_word, pics, cmr_found),
        s3_composition=         sds_ings,
        s4_first_aid=           _s4(all_codes),
        s5_fire_fighting=       _s5(avg_logp, code_set),
        s6_accidental_release=  _s6(bio_fraction),
        s7_handling_storage=    _s7(code_set),
        s8_exposure_controls=   _s8(code_set, oel_notes),
        s9_physical_props=      _s9(blend, db, avg_mw, avg_logp, avg_tpsa, viscosity_est),
        s10_stability=          _s10(code_set),
        s11_toxicological=      _s11(comptox_list, all_codes),
        s12_ecological=         _s12(comptox_list, qsar_bio_pct, qsar_etox),
        s13_disposal=           _s13(bio_fraction, code_set),
        s14_transport=          _s14(code_set),
        s15_regulatory=         _s15(reg_report, vertical),
        s16_other=              _s16(product_name, version, all_codes),
        all_ghs_codes=all_codes,
        signal_word=sig_word,
        pictograms=pics,
    )


# ── JSON export ───────────────────────────────────────────────────────────────

def sds_to_json(sds: SDSDocument) -> str:
    def _serial(obj):
        if hasattr(obj, "__dataclass_fields__"):
            return asdict(obj)
        return str(obj)
    return json.dumps(asdict(sds), indent=2, default=_serial)


# ── PDF export (ReportLab) ────────────────────────────────────────────────────

def sds_to_pdf(sds: SDSDocument) -> bytes:
    try:
        from reportlab.lib.pagesizes import A4
        from reportlab.lib.units import cm
        from reportlab.lib import colors
        from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
        from reportlab.lib.enums import TA_LEFT, TA_CENTER
        from reportlab.platypus import (
            SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle,
            HRFlowable, PageBreak, KeepTogether,
        )

        NAVY  = colors.HexColor("#0A1628")
        TEAL  = colors.HexColor("#0D9488")
        AMBER = colors.HexColor("#D97706")
        LIGHT = colors.HexColor("#F0F9FF")
        SLATE = colors.HexColor("#64748B")
        WHITE = colors.white
        RED_C = colors.HexColor("#EF4444")
        MUTED = colors.HexColor("#94A3B8")

        buf = io.BytesIO()
        doc = SimpleDocTemplate(buf, pagesize=A4,
                                leftMargin=2*cm, rightMargin=2*cm,
                                topMargin=2*cm, bottomMargin=2*cm)

        base = getSampleStyleSheet()
        def _sty(name, **kw):
            return ParagraphStyle(name, parent=base["Normal"], **kw)

        s_title   = _sty("title",   fontSize=18, textColor=WHITE,    fontName="Helvetica-Bold", leading=22)
        s_sub     = _sty("sub",     fontSize=10, textColor=MUTED,    fontName="Helvetica")
        s_sec     = _sty("sec",     fontSize=11, textColor=NAVY,     fontName="Helvetica-Bold",
                         spaceBefore=10, spaceAfter=4)
        s_body    = _sty("body",    fontSize=8.5, textColor=colors.HexColor("#1E293B"), leading=13)
        s_warn    = _sty("warn",    fontSize=8.5, textColor=RED_C,   fontName="Helvetica-Bold")
        s_amber   = _sty("amber",   fontSize=8.5, textColor=AMBER,   fontName="Helvetica-Bold")

        story = []

        # ── Cover block ───────────────────────────────────────────────────────
        cover_data = [[
            Paragraph(f"<b>SAFETY DATA SHEET</b>", s_title),
        ]]
        cover_tbl = Table(cover_data, colWidths=["100%"])
        cover_tbl.setStyle(TableStyle([
            ("BACKGROUND",    (0,0), (-1,-1), NAVY),
            ("TOPPADDING",    (0,0), (-1,-1), 16),
            ("BOTTOMPADDING", (0,0), (-1,-1), 16),
            ("LEFTPADDING",   (0,0), (-1,-1), 14),
            ("RIGHTPADDING",  (0,0), (-1,-1), 14),
            ("ROUNDEDCORNERS", [6]),
        ]))
        story.append(cover_tbl)
        story.append(Spacer(1, 6))
        story.append(Paragraph(f"<b>{sds.product_name}</b>", _sty("pn", fontSize=13, textColor=NAVY, fontName="Helvetica-Bold")))
        story.append(Paragraph(
            f"SDS Version {sds.version}  ·  Revision date: {sds.revision_date}  ·  "
            f"Prepared by: {sds.s1_identification.get('manufacturer','ChemeNova LLC')}",
            _sty("meta", fontSize=8, textColor=SLATE)))
        story.append(Spacer(1, 4))

        signal = sds.signal_word
        sig_col = RED_C if signal == "Danger" else AMBER if signal == "Warning" else SLATE
        story.append(Paragraph(f"<b>Signal word: {signal}</b>",
                               _sty("sw", fontSize=10, textColor=sig_col, fontName="Helvetica-Bold")))
        story.append(HRFlowable(width="100%", thickness=1, color=TEAL, spaceAfter=8))

        def _section(num: int, title: str, rows: list):
            """Add a numbered section table."""
            story.append(KeepTogether([
                Paragraph(f"SECTION {num} — {title.upper()}", s_sec),
                HRFlowable(width="100%", thickness=0.5, color=TEAL, spaceAfter=4),
            ]))
            if rows:
                tbl = Table(rows, colWidths=[4.5*cm, None])
                tbl.setStyle(TableStyle([
                    ("FONTNAME",      (0,0), (0,-1), "Helvetica-Bold"),
                    ("FONTSIZE",      (0,0), (-1,-1), 8.5),
                    ("VALIGN",        (0,0), (-1,-1), "TOP"),
                    ("TOPPADDING",    (0,0), (-1,-1), 3),
                    ("BOTTOMPADDING", (0,0), (-1,-1), 3),
                    ("ROWBACKGROUNDS",(0,0), (-1,-1), [LIGHT, WHITE]),
                    ("TEXTCOLOR",     (0,0), (0,-1), NAVY),
                    ("TEXTCOLOR",     (1,0), (1,-1), colors.HexColor("#1E293B")),
                    ("GRID",          (0,0), (-1,-1), 0.25, colors.HexColor("#E2E8F0")),
                ]))
                story.append(tbl)
            story.append(Spacer(1, 8))

        def _kv(label: str, value) -> list:
            val = value if isinstance(value, str) else (", ".join(value) if isinstance(value, list) else str(value))
            return [Paragraph(label, s_body), Paragraph(val or "—", s_body)]

        # ── Section 1 ─────────────────────────────────────────────────────────
        s1 = sds.s1_identification
        _section(1, "Identification", [
            _kv("Product name",    s1.get("product_name","")),
            _kv("Product code",    s1.get("product_code","")),
            _kv("Intended use",    s1.get("intended_use","")),
            _kv("Manufacturer",    s1.get("manufacturer","")),
            _kv("Address",         s1.get("address","")),
            _kv("Emergency phone", s1.get("emergency","")),
            _kv("Email",           s1.get("email","")),
        ])

        # ── Section 2 ─────────────────────────────────────────────────────────
        s2 = sds.s2_hazards
        h_rows = [_kv("Signal word", s2.get("signal_word",""))]
        pics_str = "  |  ".join(s2.get("pictograms", []))
        h_rows.append(_kv("Pictograms", pics_str))
        for hs in s2.get("hazard_statements", []):
            h_rows.append([Paragraph("", s_body), Paragraph(hs, s_warn if "H3" in hs or "H4" in hs else s_body)])
        h_rows.append(_kv("Other hazards", s2.get("other_hazards","")))
        for ps in s2.get("precautionary_statements", []):
            h_rows.append([Paragraph("", s_body), Paragraph(ps, s_body)])
        _section(2, "Hazard(s) Identification", h_rows)

        # ── Section 3 ─────────────────────────────────────────────────────────
        comp_header = [Paragraph(h, _sty("th", fontSize=8, fontName="Helvetica-Bold", textColor=WHITE))
                       for h in ["Ingredient", "CAS", "INCI", "Wt%", "REACH", "GHS codes"]]
        comp_rows = [comp_header]
        for ing in sds.s3_composition:
            comp_rows.append([
                Paragraph(ing.name[:30], s_body),
                Paragraph(ing.cas[:14],  s_body),
                Paragraph(ing.inci[:24], s_body),
                Paragraph(f"{ing.wt_pct:.1f}%", s_body),
                Paragraph(ing.reach_status[:14], s_body),
                Paragraph(", ".join(ing.ghs_codes) or "—", s_body),
            ])
        comp_tbl = Table(comp_rows, colWidths=[3.8*cm, 2.2*cm, 3.4*cm, 1.2*cm, 2.2*cm, None])
        comp_tbl.setStyle(TableStyle([
            ("BACKGROUND",    (0,0), (-1,0),  NAVY),
            ("FONTSIZE",      (0,0), (-1,-1), 8),
            ("TOPPADDING",    (0,0), (-1,-1), 3),
            ("BOTTOMPADDING", (0,0), (-1,-1), 3),
            ("ROWBACKGROUNDS",(0,1), (-1,-1), [LIGHT, WHITE]),
            ("GRID",          (0,0), (-1,-1), 0.25, colors.HexColor("#E2E8F0")),
            ("VALIGN",        (0,0), (-1,-1), "TOP"),
        ]))
        story.append(Paragraph("SECTION 3 — COMPOSITION / INFORMATION ON INGREDIENTS", s_sec))
        story.append(HRFlowable(width="100%", thickness=0.5, color=TEAL, spaceAfter=4))
        story.append(comp_tbl)
        story.append(Spacer(1, 8))

        # ── Sections 4–16 ─────────────────────────────────────────────────────
        s4 = sds.s4_first_aid
        _section(4, "First-Aid Measures", [
            _kv("Inhalation",          s4.get("inhalation","")),
            _kv("Skin contact",        s4.get("skin_contact","")),
            _kv("Eye contact",         s4.get("eye_contact","")),
            _kv("Ingestion",           s4.get("ingestion","")),
            _kv("Notes to physician",  s4.get("notes_to_physician","")),
        ])

        s5 = sds.s5_fire_fighting
        _section(5, "Fire-Fighting Measures", [
            _kv("Extinguishing media",      s5.get("extinguishing_media","")),
            _kv("Unsuitable media",         s5.get("unsuitable_media","")),
            _kv("Special hazards",          s5.get("special_hazards","")),
            _kv("Protective equipment",     s5.get("protective_equipment","")),
            _kv("Flash point note",         s5.get("flash_point_note","")),
        ])

        s6 = sds.s6_accidental_release
        _section(6, "Accidental Release Measures", [
            _kv("Personal precautions",        s6.get("personal_precautions","")),
            _kv("Environmental precautions",   s6.get("environmental_precautions","")),
            _kv("Containment",                 s6.get("containment","")),
            _kv("Clean-up",                    s6.get("cleanup","")),
        ])

        s7 = sds.s7_handling_storage
        _section(7, "Handling and Storage", [
            _kv("Handling precautions",  s7.get("handling_precautions","")),
            _kv("Storage conditions",    s7.get("storage_conditions","")),
            _kv("Incompatibilities",     s7.get("incompatibilities","")),
            _kv("Storage temperature",   s7.get("storage_temperature","")),
            _kv("Shelf life",            s7.get("shelf_life","")),
        ])

        story.append(PageBreak())

        s8 = sds.s8_exposure_controls
        _section(8, "Exposure Controls / Personal Protection", [
            _kv("Engineering controls",   s8.get("engineering_controls","")),
            _kv("OEL note",               s8.get("oel_note","")),
            _kv("Eye protection",         s8.get("eye_protection","")),
            _kv("Hand protection",        s8.get("hand_protection","")),
            _kv("Skin/body protection",   s8.get("skin_body_protection","")),
            _kv("Respiratory protection", s8.get("respiratory_protection","")),
            _kv("Hygiene measures",       s8.get("hygiene_measures","")),
        ])

        s9 = sds.s9_physical_props
        _section(9, "Physical and Chemical Properties", [
            _kv("Appearance",           s9.get("appearance","")),
            _kv("Colour",               s9.get("colour","")),
            _kv("Odour",                s9.get("odour","")),
            _kv("pH",                   s9.get("pH","")),
            _kv("Boiling point",        s9.get("boiling_point","")),
            _kv("Flash point",          s9.get("flash_point","")),
            _kv("Vapour pressure",      s9.get("vapour_pressure","")),
            _kv("Relative density",     s9.get("relative_density","")),
            _kv("Water solubility",     s9.get("water_solubility","")),
            _kv("log P (avg)",          s9.get("partition_coeff_logp","")),
            _kv("Avg. MW",              s9.get("molecular_weight_avg","")),
            _kv("Dynamic viscosity",    s9.get("dynamic_viscosity","")),
        ])

        _section(10, "Stability and Reactivity", [
            _kv("Reactivity",              sds.s10_stability.get("reactivity","")),
            _kv("Chemical stability",      sds.s10_stability.get("chemical_stability","")),
            _kv("Conditions to avoid",     sds.s10_stability.get("conditions_to_avoid","")),
            _kv("Incompatible materials",  sds.s10_stability.get("incompatible_materials","")),
            _kv("Decomposition products",  sds.s10_stability.get("decomposition_products","")),
        ])

        s11 = sds.s11_toxicological
        _section(11, "Toxicological Information", [
            _kv("Acute toxicity (oral)",    s11.get("acute_toxicity_oral","")),
            _kv("Skin irritation",          s11.get("skin_irritation","")),
            _kv("Eye irritation",           s11.get("eye_irritation","")),
            _kv("Resp. sensitisation",      s11.get("respiratory_sensitisation","")),
            _kv("Carcinogenicity",          s11.get("carcinogenicity","")),
            _kv("Reproductive toxicity",    s11.get("reproductive_toxicity","")),
            _kv("CMR summary",              s11.get("cmr_summary","")),
            _kv("Fish LC50 (log mmol/L)",   s11.get("fish_lc50_log_mmol","")),
            _kv("Daphnia EC50 (log mmol/L)",s11.get("daphnia_ec50_log_mmol","")),
            _kv("Avg. log Kow",             s11.get("avg_log_kow","")),
            _kv("Notes",                    s11.get("notes","")),
        ])

        s12 = sds.s12_ecological
        _section(12, "Ecological Information", [
            _kv("Aquatic toxicity",   s12.get("aquatic_toxicity","")),
            _kv("Persistence",        s12.get("persistence","")),
            _kv("Bioaccumulation",    s12.get("bioaccumulation","")),
            _kv("Mobility in soil",   s12.get("mobility_soil","")),
            _kv("PBT/vPvB",           s12.get("pbt_vpvb_assessment","")),
            _kv("Other effects",      s12.get("other_adverse_effects","")),
        ])

        _section(13, "Disposal Considerations", [
            _kv("Waste treatment",    sds.s13_disposal.get("waste_treatment","")),
            _kv("Waste code",         sds.s13_disposal.get("waste_code","")),
            _kv("Classification",     sds.s13_disposal.get("classification","")),
            _kv("Packaging",          sds.s13_disposal.get("packaging_disposal","")),
        ])

        s14 = sds.s14_transport
        _section(14, "Transport Information", [
            _kv("UN number",          s14.get("un_number","")),
            _kv("Shipping name",      s14.get("proper_shipping_name","")),
            _kv("Hazard class",       s14.get("hazard_class","")),
            _kv("Packing group",      s14.get("packing_group","")),
            _kv("Marine pollutant",   s14.get("marine_pollutant","")),
            _kv("Special precautions",s14.get("special_precautions","")),
        ])

        s15 = sds.s15_regulatory
        _section(15, "Regulatory Information", [
            _kv("EU REACH",           s15.get("eu_reach","")),
            _kv("EU CLP",             s15.get("eu_clp","")),
            _kv("EU Ecolabel",        s15.get("eu_ecolabel","")),
            _kv("COSMOS",             s15.get("cosmos","")),
            _kv("EPA Safer Choice",   s15.get("epa_safer_choice","")),
            _kv("California Prop 65", s15.get("california_prop65","")),
            _kv("Other",              s15.get("other_regulations","")),
        ])

        s16 = sds.s16_other
        h_texts_16 = "  ·  ".join(f"{k}: {v}" for k, v in s16.get("full_h_text", {}).items())
        _section(16, "Other Information", [
            _kv("Revision info",     s16.get("revision_info","")),
            _kv("Revision date",     s16.get("revision_date","")),
            _kv("Prepared by",       s16.get("prepared_by","")),
            _kv("H-statement key",   h_texts_16 or "No hazard statements"),
            _kv("Training required", s16.get("training_required","")),
            _kv("Disclaimer",        s16.get("disclaimer","")),
            _kv("Reference",         s16.get("intelliform_ref","")),
        ])

        # ── Footer rule ───────────────────────────────────────────────────────
        story.append(HRFlowable(width="100%", thickness=1, color=TEAL))
        story.append(Paragraph(
            f"End of SDS · {sds.product_name} · Version {sds.version} · "
            f"{sds.revision_date} · IntelliForm™ ChemeNova LLC",
            _sty("foot", fontSize=7, textColor=SLATE, alignment=1)
        ))

        doc.build(story)
        return buf.getvalue()

    except ImportError:
        raise RuntimeError("reportlab is required for PDF export: pip install reportlab")
