"""
modules/pharma.py
Pharmaceutical Deep Dive Module for IntelliForm v1.1

Implements:
  1. BCS Classification (Biopharmaceutics Classification System)
  2. API-Excipient Compatibility Matrix
  3. ICH Q1A Stability Zone Predictions
  4. Dosage Form Engineering (tablet, capsule, liquid, topical)
  5. Excipient Function Optimization
  6. Regulatory Pathway Guidance (NDA/ANDA/505(b)(2))
  7. Manufacturing Route Intelligence
  8. Drug-Excipient Interaction Flags

References:
  - ICH Q8(R2) Pharmaceutical Development
  - ICH Q1A(R2) Stability Testing
  - BCS Classification: Amidon et al., Pharm Res 1995
  - FDA GRASE list / USP-NF monographs
  - EMA Excipient Database
"""
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple
import pandas as pd


# ── BCS Classification ─────────────────────────────────────────────────────────

@dataclass
class BCSProfile:
    bcs_class: str          # I, II, III, IV
    solubility: str         # High / Low
    permeability: str       # High / Low
    absorption_prediction: str
    formulation_strategy: str
    enabling_excipients: List[str]
    bioavailability_risk: str  # Low / Medium / High
    ivivc_potential: str       # Good / Moderate / Poor


BCS_STRATEGIES = {
    "I": BCSProfile(
        bcs_class="I",
        solubility="High", permeability="High",
        absorption_prediction="Good absorption expected — dissolution rarely limiting",
        formulation_strategy="Immediate release formulation. Simple direct compression viable. "
                             "Focus on achieving rapid dissolution >85% in 15 min.",
        enabling_excipients=["Microcrystalline Cellulose", "Lactose Monohydrate",
                             "Croscarmellose Sodium", "Magnesium Stearate"],
        bioavailability_risk="Low",
        ivivc_potential="Good — Level A IVIVC achievable",
    ),
    "II": BCSProfile(
        bcs_class="II",
        solubility="Low", permeability="High",
        absorption_prediction="Absorption limited by solubility/dissolution rate — formulation critical",
        formulation_strategy="Solubility/dissolution enhancement required. Consider: "
                             "(1) Solid dispersion/HME with PVP or HPMC-AS, "
                             "(2) Nanoparticle/nanosuspension, "
                             "(3) Lipid-based systems (SMEDDS/SEDDS), "
                             "(4) Cyclodextrin complexation, "
                             "(5) Salt formation if ionizable.",
        enabling_excipients=["Povidone (PVP K30)", "HPMC-AS", "HP-Beta-Cyclodextrin",
                             "Poloxamer 407", "Cremophor EL", "Medium Chain Triglycerides (pharma)",
                             "Tween 80 (Polysorbate 80)", "Vitamin E TPGS"],
        bioavailability_risk="High — formulation-dependent",
        ivivc_potential="Good if dissolution-limited — Level A IVIVC possible",
    ),
    "III": BCSProfile(
        bcs_class="III",
        solubility="High", permeability="Low",
        absorption_prediction="Dissolution rapid but permeability limits absorption",
        formulation_strategy="Permeation enhancement strategies: "
                             "(1) Tight junction modulators (chitosan, EDTA), "
                             "(2) Efflux pump inhibitors (Tween 80, Vitamin E TPGS), "
                             "(3) Mucoadhesive systems for extended contact time, "
                             "(4) Nanocarrier encapsulation.",
        enabling_excipients=["Chitosan (pharma grade)", "Polycarbophil",
                             "Carbopol 934 (pharma)", "Tween 80 (Polysorbate 80)",
                             "Sodium EDTA", "Hypromellose (HPMC E5)"],
        bioavailability_risk="Medium — permeability dependent",
        ivivc_potential="Poor — IVIVC difficult due to permeability limitation",
    ),
    "IV": BCSProfile(
        bcs_class="IV",
        solubility="Low", permeability="Low",
        absorption_prediction="Both solubility and permeability limit absorption — challenging",
        formulation_strategy="Most challenging class. Combined approach required: "
                             "(1) Lipid-based self-emulsifying systems (SEDDS/SMEDDS) — best approach, "
                             "(2) Nanoparticle + absorption enhancer combination, "
                             "(3) Prodrug strategy if feasible, "
                             "(4) Parenteral/alternative route consideration.",
        enabling_excipients=["Medium Chain Triglycerides (pharma)", "Cremophor EL",
                             "Tween 80 (Polysorbate 80)", "Lecithin (IV grade)",
                             "Poloxamer 407", "Ethanol (pharma)"],
        bioavailability_risk="High — highly formulation-dependent",
        ivivc_potential="Poor — route-dependent behavior",
    ),
}


# ── API-Excipient Compatibility ────────────────────────────────────────────────

@dataclass
class CompatibilityResult:
    ingredient: str
    compatibility_class: str   # Compatible / Conditional / Incompatible / Unknown
    interaction_type: str      # None / Chemical / Physical / Adsorption / Catalytic
    mechanism: str
    severity: str              # None / Mild / Moderate / Severe
    mitigation: str
    literature_ref: str


# Known incompatibility patterns — function-based rules
INCOMPATIBILITY_RULES = [
    {
        "trigger_functions": ["Pharma Lubricant"],
        "trigger_ingredients": ["Magnesium Stearate"],
        "condition": "hydrophobic lubricant",
        "interaction": "Chemical/Physical",
        "mechanism": "Magnesium stearate hydrophobic film reduces dissolution rate. "
                    "Over-lubrication retards disintegration. Interacts with amine APIs.",
        "severity": "Moderate",
        "mitigation": "Limit to ≤0.5% for immediate release. Use alternative lubricants "
                     "(sodium stearyl fumarate, glyceryl behenate) for amine APIs.",
        "ref": "Bolhuis & Chowhan, Tablets and Capsules, 1996",
    },
    {
        "trigger_functions": ["Pharma Sweetener"],
        "trigger_ingredients": ["Lactose Monohydrate"],
        "condition": "reducing sugar",
        "interaction": "Chemical",
        "mechanism": "Maillard reaction with primary amine APIs. Browning and potency loss. "
                    "Accelerated at elevated temperature/humidity.",
        "severity": "Severe",
        "mitigation": "Avoid lactose with amine-containing APIs. Use mannitol, "
                     "microcrystalline cellulose, or calcium phosphate as alternatives.",
        "ref": "Wirth et al., Drug Dev Ind Pharm, 1998",
    },
    {
        "trigger_functions": ["Pharma Disintegrant"],
        "trigger_ingredients": ["Crospovidone", "Croscarmellose Sodium"],
        "condition": "high-dose disintegrant",
        "interaction": "Adsorption",
        "mechanism": "Crosslinked polymers may adsorb basic APIs reducing dissolution. "
                    "More pronounced with cationic drugs at physiological pH.",
        "severity": "Mild",
        "mitigation": "Evaluate API adsorption in dissolution media. Consider external "
                     "addition of disintegrant. Test at 2%, 4%, 6% levels.",
        "ref": "Kornblum & Stoopak, J Pharm Sci, 1973",
    },
    {
        "trigger_functions": ["Pharma Antacid"],
        "trigger_ingredients": ["Magnesium Oxide (pharma)", "Calcium Carbonate (pharma)",
                                "Aluminum Hydroxide (pharma)"],
        "condition": "alkaline excipient",
        "interaction": "Chemical",
        "mechanism": "Alkaline conditions can catalyze hydrolysis of ester and amide APIs. "
                    "pH elevation may reduce solubility of basic APIs.",
        "severity": "Moderate",
        "mitigation": "Monitor pH of blend. Use buffering agents. "
                     "Evaluate stability at formulation pH.",
        "ref": "Connors et al., Chemical Stability of Pharmaceuticals, 1986",
    },
    {
        "trigger_functions": ["Pharma Flow Aid", "Pharma Glidant"],
        "trigger_ingredients": ["Silicon Dioxide (Aerosil)", "Colloidal Silicon Dioxide"],
        "condition": "high surface area silica",
        "interaction": "Adsorption",
        "mechanism": "High surface area silica can adsorb drug molecules reducing bioavailability. "
                    "Particularly problematic with BCS Class II/IV APIs.",
        "severity": "Mild",
        "mitigation": "Limit to ≤0.5%. Consider hydrophobic grade (R972) "
                     "for moisture-sensitive APIs.",
        "ref": "Rao et al., Int J Pharm, 2001",
    },
    {
        "trigger_functions": ["Pharma Enteric Polymer", "Pharma Sustained Release"],
        "trigger_ingredients": ["Eudragit L100", "Eudragit RS"],
        "condition": "acrylic polymer",
        "interaction": "Physical",
        "mechanism": "Acrylic polymers require plasticizer for proper film formation. "
                    "Without plasticizer: brittle coat, cracking, dose dumping risk.",
        "severity": "Moderate",
        "mitigation": "Include triethyl citrate (15-20% of polymer) or TEC as plasticizer. "
                     "Cure coating at 40°C for 24h to complete coalescence.",
        "ref": "Lehmann, Drugs Made in Germany, 1989",
    },
    {
        "trigger_functions": ["Pharma Preservative"],
        "trigger_ingredients": ["Benzalkonium Chloride"],
        "condition": "cationic preservative",
        "interaction": "Chemical",
        "mechanism": "Benzalkonium chloride incompatible with anionic APIs and excipients. "
                    "Forms precipitates with soaps, anionic surfactants, citrates, nitrates.",
        "severity": "Severe",
        "mitigation": "Avoid with anionic drugs. Use phenoxyethanol or sorbate-based "
                     "preservative systems for anionic API formulations.",
        "ref": "Lam et al., J Pharm Sci, 2012",
    },
    {
        "trigger_functions": ["Pharma Solubilizer", "Pharma Solubilizer/Complexer"],
        "trigger_ingredients": ["HP-Beta-Cyclodextrin", "Hydroxypropyl Cellulose (HPC)"],
        "condition": "complexation agent",
        "interaction": "Physical",
        "mechanism": "Cyclodextrins form inclusion complexes — can compete with other "
                    "excipients for API. May alter dissolution rate and absorption.",
        "severity": "Mild",
        "mitigation": "Characterize complex stoichiometry. Phase solubility studies "
                     "required. Monitor for competitor displacement.",
        "ref": "Loftsson & Brewster, J Pharm Pharmacol, 1996",
    },
]


# ── ICH Stability Zones ────────────────────────────────────────────────────────

@dataclass
class ICHStabilityZone:
    zone: str
    regions: List[str]
    long_term: str
    intermediate: str
    accelerated: str
    critical_conditions: str
    packaging_recommendation: str


ICH_STABILITY_ZONES = {
    "I": ICHStabilityZone(
        zone="I — Temperate",
        regions=["USA", "Canada", "Northern Europe", "Japan"],
        long_term="25°C/60% RH for 12 months minimum",
        intermediate="Not required",
        accelerated="40°C/75% RH for 6 months",
        critical_conditions="Low humidity stress — hydrolysis less likely",
        packaging_recommendation="HDPE or glass with desiccant for moisture-sensitive APIs",
    ),
    "II": ICHStabilityZone(
        zone="II — Mediterranean/Subtropical",
        regions=["Southern Europe", "Japan", "USA (some regions)"],
        long_term="25°C/60% RH for 12 months",
        intermediate="30°C/65% RH for 6 months",
        accelerated="40°C/75% RH for 6 months",
        critical_conditions="Intermediate humidity — photodegradation risk in summer",
        packaging_recommendation="Light-protective packaging, desiccant for hygroscopic formulations",
    ),
    "III": ICHStabilityZone(
        zone="III — Hot/Dry",
        regions=["Middle East", "North Africa", "parts of Asia"],
        long_term="30°C/35% RH for 12 months",
        intermediate="Not required",
        accelerated="40°C/NMT 25% RH for 6 months",
        critical_conditions="High temperature stress — oxidation and thermal degradation",
        packaging_recommendation="Hermetically sealed packaging, antioxidants in formulation",
    ),
    "IVA": ICHStabilityZone(
        zone="IVA — Hot/Humid",
        regions=["Sub-Saharan Africa", "Southeast Asia", "India", "South America"],
        long_term="30°C/65% RH for 12 months",
        intermediate="Not required",
        accelerated="40°C/75% RH for 6 months",
        critical_conditions="High temperature + humidity — hydrolysis, microbial growth",
        packaging_recommendation="Moisture-barrier packaging (Alu-Alu blisters), "
                                "antimicrobial preservatives for liquids",
    ),
    "IVB": ICHStabilityZone(
        zone="IVB — Hot/Very Humid",
        regions=["Brazil", "parts of SE Asia", "tropical regions"],
        long_term="30°C/75% RH for 12 months",
        intermediate="Not required",
        accelerated="40°C/75% RH for 6 months",
        critical_conditions="Maximum humidity stress — aggressive hydrolysis and microbial risk",
        packaging_recommendation="Tropical packaging — Alu-Alu blisters or HDPE with strong desiccant",
    ),
}


# ── Dosage Form Profiles ───────────────────────────────────────────────────────

@dataclass
class DosageFormProfile:
    form: str
    description: str
    typical_composition: Dict[str, Tuple[float, float]]   # function: (min%, max%)
    key_excipients: List[str]
    critical_quality_attributes: List[str]
    manufacturing_process: str
    regulatory_considerations: List[str]
    bcs_suitability: List[str]   # Which BCS classes this form suits


DOSAGE_FORMS = {
    "immediate_release_tablet": DosageFormProfile(
        form="Immediate Release Tablet",
        description="Conventional oral solid dosage form for systemic delivery. "
                   "Most common, well-understood, cost-effective.",
        typical_composition={
            "Filler/Diluent": (50.0, 80.0),
            "Binder": (2.0, 10.0),
            "Disintegrant": (2.0, 8.0),
            "Lubricant": (0.25, 1.5),
            "Glidant": (0.1, 0.5),
            "API": (0.5, 30.0),
        },
        key_excipients=[
            "Microcrystalline Cellulose PH102 (filler/binder)",
            "Lactose Monohydrate (filler — avoid with amine APIs)",
            "Croscarmellose Sodium (disintegrant)",
            "Magnesium Stearate (lubricant)",
            "Colloidal Silicon Dioxide (glidant)",
            "Povidone (PVP K30) (wet granulation binder)",
        ],
        critical_quality_attributes=[
            "Disintegration time (<15 min for IR)",
            "Dissolution (>85% in 30 min for BCS Class I)",
            "Hardness (4-20 kP)",
            "Friability (<1%)",
            "Weight uniformity (±5%)",
            "Content uniformity (RSD <2%)",
        ],
        manufacturing_process="Direct Compression (preferred for MCC-dominant blends) or "
                             "Wet Granulation (for poor-flow APIs) or "
                             "Dry Granulation/Roller Compaction (for moisture-sensitive APIs)",
        regulatory_considerations=[
            "ICH Q8 — Design Space definition",
            "ICH Q9 — Quality Risk Management",
            "ICH Q10 — Pharmaceutical Quality System",
            "FDA Q1 — Dissolution method development",
            "21 CFR 211 — cGMP compliance",
        ],
        bcs_suitability=["I", "III"],
    ),

    "modified_release_tablet": DosageFormProfile(
        form="Modified Release Tablet (SR/CR/ER)",
        description="Extended, sustained, or controlled release tablet. "
                   "Reduces dosing frequency, maintains therapeutic levels.",
        typical_composition={
            "Matrix Polymer": (10.0, 40.0),
            "Filler": (20.0, 60.0),
            "Lubricant": (0.25, 1.5),
            "Film Coat": (2.0, 5.0),
            "API": (1.0, 40.0),
        },
        key_excipients=[
            "Methocel E5 / K100M (HPMC matrix — hydrophilic)",
            "Eudragit RS (acrylic matrix — pH-independent)",
            "Eudragit L100 (enteric coat — pH >6)",
            "Ethylcellulose (hydrophobic matrix)",
            "Glyceryl Behenate (lipid matrix)",
            "Triethyl Citrate (plasticizer for Eudragit)",
        ],
        critical_quality_attributes=[
            "Dissolution profile at multiple pH (1.2, 4.5, 6.8)",
            "Drug release mechanism (zero-order vs first-order)",
            "Food effect evaluation",
            "Dose dumping assessment",
            "Polymer hydration rate (for HPMC matrix)",
        ],
        manufacturing_process="Matrix Tablet: Direct compression or granulation with polymer. "
                             "Coated Pellets: Fluid-bed coating with functional polymer. "
                             "Osmotic: Complex manufacturing, high cost.",
        regulatory_considerations=[
            "FDA Guidance for Industry: Modified Release Solid Oral Dosage Forms",
            "IVIVC required for Level A — NDA benefit",
            "Biowaivers not available for MR formulations",
            "ICH Q8 — Enhanced approach strongly recommended",
        ],
        bcs_suitability=["I", "II"],
    ),

    "hard_capsule": DosageFormProfile(
        form="Hard Capsule (HPMC or Gelatin)",
        description="Powder or pellet fill in hard shell capsule. "
                   "Easier formulation development than tablets, no compression needed.",
        typical_composition={
            "API + Filler blend": (80.0, 95.0),
            "Capsule Shell": (5.0, 20.0),
            "Lubricant": (0.25, 1.0),
            "Disintegrant (if needed)": (0.0, 5.0),
        },
        key_excipients=[
            "HPMC Capsule Shell (veg — suitable for halal/kosher)",
            "Gelatin Capsule Shell (animal-derived)",
            "Microcrystalline Cellulose (filler)",
            "Mannitol DC Grade (filler — low hygroscopicity)",
            "Magnesium Stearate (lubricant — max 0.5%)",
            "Silicon Dioxide (glidant)",
        ],
        critical_quality_attributes=[
            "Capsule fill weight uniformity",
            "Disintegration (<30 min for gelatin, <60 min for HPMC)",
            "Moisture content of fill (critical for gelatin shells)",
            "Shell integrity (brittleness at low RH)",
        ],
        manufacturing_process="Fill powder blend into capsule shells. "
                             "Pellet fill: fluid-bed granulation + coating + filling. "
                             "Liquid-fill: hot-melt fill with lipid systems.",
        regulatory_considerations=[
            "HPMC shell: no bovine/porcine concerns — global market access",
            "Gelatin shell: BSE/TSE certification required for EU",
            "Shell coloring: EU approved colorants only",
            "Shell size selection based on fill density and dose",
        ],
        bcs_suitability=["I", "II", "III"],
    ),

    "oral_solution": DosageFormProfile(
        form="Oral Solution / Suspension",
        description="Liquid dosage form for pediatrics, dysphagia, or high-dose APIs. "
                   "Bioavailability advantages for poorly soluble drugs.",
        typical_composition={
            "Solvent/Cosolvent": (20.0, 70.0),
            "Sweetener/Taste Masking": (5.0, 30.0),
            "Preservative": (0.05, 0.5),
            "Buffer/pH Modifier": (0.1, 2.0),
            "Viscosity Modifier": (0.1, 3.0),
            "API": (0.1, 20.0),
        },
        key_excipients=[
            "Sorbitol Solution 70% (sweetener/vehicle)",
            "Glycerin (USP grade) (cosolvent/sweetener)",
            "Propylene Glycol (USP) (cosolvent)",
            "Sodium Benzoate (preservative)",
            "Citric Acid Monohydrate (buffer/flavor)",
            "Sucrose (pharma) (sweetener)",
            "Xanthan Gum (viscosity modifier for suspension)",
        ],
        critical_quality_attributes=[
            "pH (typically 3-7 for stability)",
            "Osmolality (iso-osmotic preferred for pediatrics)",
            "Preservative efficacy (USP <51>)",
            "Chemical stability of API in solution",
            "Microbial limits (USP <61>/<62>)",
            "Particle size for suspensions (D90 < 10 µm)",
        ],
        manufacturing_process="Dissolve in order of decreasing polarity. "
                             "Buffer before API. Preserve. QS to volume. "
                             "Aseptic processing if sterile required.",
        regulatory_considerations=[
            "Pediatric dose verification (FDA Pediatric Rule)",
            "Alcohol content limits for pediatric products",
            "Preservative compatibility with container",
            "Reconstitution instructions for dry powder",
        ],
        bcs_suitability=["I", "II", "III", "IV"],
    ),

    "topical_semisolid": DosageFormProfile(
        form="Topical Semisolid (Cream/Gel/Ointment)",
        description="Topical application for local or transdermal drug delivery.",
        typical_composition={
            "Emollient/Oil Phase": (10.0, 40.0),
            "Emulsifier": (2.0, 10.0),
            "Thickener/Gelling Agent": (0.5, 5.0),
            "Preservative": (0.05, 0.5),
            "Water/Aqueous Phase": (30.0, 70.0),
            "API": (0.1, 10.0),
        },
        key_excipients=[
            "White Petrolatum (pharma) (ointment base)",
            "Cetearyl Alcohol (emollient/emulsifier)",
            "Ceteareth-20 (o/w emulsifier)",
            "Carbopol 934 (pharma) (gel base)",
            "Methylparaben + Propylparaben (preservative system)",
            "Propylene Glycol (USP) (humectant/penetration enhancer)",
            "Triethanolamine (pH neutralizer for Carbopol)",
        ],
        critical_quality_attributes=[
            "Viscosity/rheology (spreadability)",
            "pH (4.5-7.5 for skin compatibility)",
            "Preservative efficacy (USP <51>)",
            "Particle size for suspensions",
            "In vitro drug release (Franz diffusion cell)",
            "Physical stability (freeze-thaw cycling)",
        ],
        manufacturing_process="Cream: heat aqueous and oil phases separately to 70-75°C, "
                             "add aqueous to oil phase under homogenization, cool with mixing. "
                             "Gel: disperse polymer cold, neutralize, add API.",
        regulatory_considerations=[
            "FDA Q2: Analytical procedure validation",
            "ICH Q3B: Degradation products",
            "Semi-solid FDA guidance on ANDA requirements",
            "Preservative efficacy testing mandatory",
        ],
        bcs_suitability=["I", "II", "III", "IV"],
    ),

    "parenteral": DosageFormProfile(
        form="Parenteral (IV/IM/SC Injection)",
        description="Sterile injectable product for IV, IM, or SC administration. "
                   "100% bioavailability. Most regulated dosage form.",
        typical_composition={
            "Water for Injection": (80.0, 99.0),
            "Tonicity Agent": (0.5, 2.0),
            "Buffer": (0.1, 1.0),
            "Solubilizer/Surfactant": (0.0, 5.0),
            "Preservative (MDV only)": (0.0, 0.5),
            "API": (0.1, 15.0),
        },
        key_excipients=[
            "Water for Injection (WFI) — sterile, pyrogen-free",
            "Sodium Chloride (tonicity)",
            "Sodium Phosphate Dibasic + Citric Acid (buffer)",
            "Polysorbate 80 (solubilizer)",
            "Benzyl Alcohol (preservative — MDV only)",
            "Mannitol (isotonicity + lyophilization bulking agent)",
            "Lecithin (IV grade) (emulsifier for lipid emulsions)",
        ],
        critical_quality_attributes=[
            "Sterility (USP <71>)",
            "Bacterial endotoxins (USP <85>) — <5 EU/kg",
            "Particulate matter (USP <788>/<789>)",
            "pH (4-9 for IV; 3-10 for IM/SC)",
            "Osmolality (270-310 mOsm/kg for IV)",
            "Extractables and leachables from container",
        ],
        manufacturing_process="Aseptic processing or terminal sterilization (preferred). "
                             "Terminal: autoclave 121°C/15min for heat-stable APIs. "
                             "Aseptic: sterile filtration (0.22 µm) + fill in ISO 5. "
                             "Lyophilization for unstable APIs in aqueous solution.",
        regulatory_considerations=[
            "21 CFR 211 Subpart K — Sterile drug products",
            "ICH Q5C — Stability of biotechnology products",
            "FDA Guidance: Sterile Drug Products (aseptic processing)",
            "Container closure system integrity testing required",
            "Extractables/leachables per USP <661>",
        ],
        bcs_suitability=["I", "II", "III", "IV"],
    ),
}


# ── Excipient Function Optimizer ──────────────────────────────────────────────

@dataclass
class ExcipientRecommendation:
    function: str
    primary_choice: str
    alternative_1: str
    alternative_2: str
    rationale: str
    typical_loading: str
    compatibility_note: str


EXCIPIENT_RECOMMENDATIONS = {
    "filler_general": ExcipientRecommendation(
        function="Filler/Diluent — General",
        primary_choice="Microcrystalline Cellulose PH102",
        alternative_1="Mannitol DC Grade",
        alternative_2="Calcium Phosphate Dibasic",
        rationale="MCC offers excellent compressibility, disintegration contribution, "
                 "and chemical inertness. Most versatile pharmaceutical filler.",
        typical_loading="40-80% of tablet weight",
        compatibility_note="Avoid lactose with amine-containing APIs (Maillard reaction)",
    ),
    "filler_moisture_sensitive": ExcipientRecommendation(
        function="Filler — Moisture-Sensitive API",
        primary_choice="Mannitol DC Grade",
        alternative_1="Calcium Phosphate Dibasic",
        alternative_2="Calcium Carbonate (pharma)",
        rationale="Mannitol: low moisture content, no reducing sugar, good flow. "
                 "Ideal for moisture-sensitive and amine APIs.",
        typical_loading="40-70% of tablet weight",
        compatibility_note="Mannitol may absorb moisture if not properly sealed — "
                          "use with desiccant packaging",
    ),
    "disintegrant": ExcipientRecommendation(
        function="Disintegrant",
        primary_choice="Croscarmellose Sodium",
        alternative_1="Sodium Starch Glycolate (type A)",
        alternative_2="Crospovidone",
        rationale="Croscarmellose: best swelling capacity, rapid disintegration. "
                 "SSG: excellent for starch-based systems. Crospovidone: wicking mechanism.",
        typical_loading="2-8% (intragranular) or 2-4% (extragranular)",
        compatibility_note="Crospovidone may adsorb basic APIs — "
                          "evaluate if BCS II/IV drug",
    ),
    "lubricant_standard": ExcipientRecommendation(
        function="Lubricant — Standard",
        primary_choice="Magnesium Stearate",
        alternative_1="Sodium Stearyl Fumarate",
        alternative_2="Glyceryl Behenate",
        rationale="Magnesium stearate: most cost-effective, widely used. "
                 "Sodium stearyl fumarate: lower sensitivity to over-lubrication. "
                 "Glyceryl behenate: matrix-forming lubricant for SR.",
        typical_loading="0.25-1.5% (MgSt); 1-2% (SSF); 1-5% (GB)",
        compatibility_note="MgSt + amine APIs = incompatibility risk. "
                          "Use SSF or GB for amine-containing drugs.",
    ),
    "binder_wet_granulation": ExcipientRecommendation(
        function="Binder — Wet Granulation",
        primary_choice="Povidone (PVP K30)",
        alternative_1="Hypromellose (HPMC E5)",
        alternative_2="Hydroxypropyl Cellulose (HPC)",
        rationale="PVP K30: excellent binding, soluble in water and alcohols, "
                 "compatible with most APIs. HPMC E5: natural-origin, suitable for "
                 "film coating and granulation. HPC: good thermoplastic properties.",
        typical_loading="2-8% (PVP); 2-6% (HPMC); 2-5% (HPC)",
        compatibility_note="PVP is hygroscopic — ensure proper storage and "
                          "moisture content control in granules",
    ),
    "film_coat_ir": ExcipientRecommendation(
        function="Film Coat — Immediate Release",
        primary_choice="Opadry (HPMC film coat)",
        alternative_1="Hypromellose (HPMC E5)",
        alternative_2="Kollicoat IR",
        rationale="Opadry: complete coating system, color + HPMC + plasticizer. "
                 "Plain HPMC: flexible, low-cost, broad compatibility. "
                 "Kollicoat: PVA-based, rapid dissolution, smooth finish.",
        typical_loading="2-5% weight gain",
        compatibility_note="Aqueous film coating preferred. "
                          "Add 10-15% plasticizer (PEG 400 or triethyl citrate) "
                          "for plain HPMC systems.",
    ),
    "enteric_coat": ExcipientRecommendation(
        function="Enteric Coating",
        primary_choice="Eudragit L100",
        alternative_1="Hypromellose Acetate Succinate (HPMC-AS)",
        alternative_2="Shellac (pharmaceutical)",
        rationale="Eudragit L100: dissolves at pH >6.0, robust coating, "
                 "well-characterized. HPMC-AS: natural polymer, oral bioavailability "
                 "enhancer for BCS II. Shellac: natural, dissolves >7.",
        typical_loading="8-12% weight gain",
        compatibility_note="Eudragit requires triethyl citrate as plasticizer (15-20% of polymer). "
                          "Cure at 40°C/24h mandatory to prevent dose dumping.",
    ),
    "preservative_liquid": ExcipientRecommendation(
        function="Preservative — Oral Liquid",
        primary_choice="Sodium Benzoate",
        alternative_1="Potassium Sorbate",
        alternative_2="Methyl Paraben (pharma)",
        rationale="Sodium benzoate: effective at pH <4.5, GRAS, cost-effective. "
                 "Potassium sorbate: effective pH 3-6, broader spectrum. "
                 "Methylparaben: broad spectrum, requires propylparaben synergy.",
        typical_loading="0.1-0.5% (benzoate); 0.1-0.2% (sorbate); "
                       "0.05-0.2% methyl + 0.01-0.02% propyl paraben",
        compatibility_note="Benzoate + acidic APIs may cause precipitation. "
                          "Paraben adsorption to polyols (sorbitol, PG) reduces efficacy — "
                          "increase loading or use alternative.",
    ),
}


# ── Regulatory Pathway Guide ──────────────────────────────────────────────────

@dataclass
class RegulatoryPathway:
    pathway: str
    description: str
    timeline: str
    cost_estimate: str
    key_studies: List[str]
    advantages: List[str]
    requirements: List[str]


REGULATORY_PATHWAYS = {
    "NDA_505b1": RegulatoryPathway(
        pathway="NDA 505(b)(1) — Full New Drug Application",
        description="Complete safety and efficacy data package for new chemical entity (NCE). "
                   "Required for novel APIs not previously approved.",
        timeline="10-15 years (Phase I-III + FDA review)",
        cost_estimate="$1.5B - $2.5B average",
        key_studies=["Phase I (safety/PK)", "Phase II (efficacy signal)",
                    "Phase III (pivotal efficacy + safety)", "Long-term stability (ICH Q1A)",
                    "Carcinogenicity studies", "Reproductive toxicology"],
        advantages=["Full IP protection (5-year NCE exclusivity)",
                   "Brand pricing power", "First-mover advantage"],
        requirements=["Complete preclinical package", "Clinical trial data",
                     "Chemistry/Manufacturing/Controls (CMC)", "Risk Evaluation and Mitigation Strategy (REMS)"],
    ),
    "NDA_505b2": RegulatoryPathway(
        pathway="NDA 505(b)(2) — Hybrid Application",
        description="Relies partially on published literature or FDA's finding of safety "
                   "and efficacy for a previously approved drug. New formulation, new route, "
                   "new indication of approved API.",
        timeline="2-5 years",
        cost_estimate="$10M - $100M",
        key_studies=["Bioequivalence or comparative BA studies",
                    "Bridging studies", "Stability per ICH Q1A",
                    "Formulation characterization (ICH Q8)"],
        advantages=["Faster timeline vs full NDA", "Lower cost", "3-year exclusivity possible",
                   "Pediatric exclusivity (6 months) if pediatric studies done"],
        requirements=["Demonstration of reliance on existing data",
                     "Patent certification (Paragraph I-IV)", "Bioequivalence data"],
    ),
    "ANDA_505j": RegulatoryPathway(
        pathway="ANDA 505(j) — Abbreviated New Drug Application (Generic)",
        description="Generic drug application demonstrating bioequivalence to Reference Listed Drug (RLD). "
                   "No new clinical efficacy data required.",
        timeline="2-4 years",
        cost_estimate="$1M - $5M",
        key_studies=["Bioequivalence study (fed/fasted)", "Dissolution method development",
                    "Stability (ICH Q1A) — accelerated 6 months + ongoing",
                    "Impurity profiling per ICH Q3B",
                    "Container closure integrity"],
        advantages=["Lowest cost regulatory pathway", "Immediate generic market entry",
                   "FDA priority review incentives available"],
        requirements=["Pharmaceutical equivalence to RLD",
                     "Bioequivalence (90% CI 80-125%)",
                     "Same dosage form, strength, route",
                     "Patent clearance"],
    ),
}


# ── Main Pharma Deep Dive Engine ──────────────────────────────────────────────

@dataclass
class PharmaDeepDiveResult:
    # BCS
    bcs_class: Optional[str]
    bcs_profile: Optional[BCSProfile]

    # Dosage form
    recommended_dosage_form: str
    dosage_form_profile: Optional[DosageFormProfile]
    alternative_forms: List[str]

    # Excipient recommendations
    excipient_recommendations: List[ExcipientRecommendation]

    # Compatibility
    compatibility_results: List[CompatibilityResult]
    compatibility_summary: str

    # ICH stability
    stability_zone: str
    stability_profile: Optional[ICHStabilityZone]
    stability_concerns: List[str]
    packaging_recommendation: str

    # Manufacturing
    manufacturing_route: str
    manufacturing_rationale: str

    # Regulatory pathway
    suggested_pathway: str
    pathway_profile: Optional[RegulatoryPathway]

    # Overall
    development_risks: List[str]
    development_recommendations: List[str]
    icm_score: float  # IntelliForm Complexity Metric 0-100


def run_pharma_deep_dive(
    blend: Dict[str, float],
    db: pd.DataFrame,
    bcs_class: str = "I",
    dosage_form: str = "immediate_release_tablet",
    target_markets: List[str] = None,
    is_generic: bool = False,
    is_pediatric: bool = False,
) -> PharmaDeepDiveResult:
    """
    Run full pharmaceutical deep dive analysis on a blend.
    """
    if target_markets is None:
        target_markets = ["USA", "EU"]

    idx = db.set_index("Ingredient") if "Ingredient" in db.columns else db

    # BCS Profile
    bcs_profile = BCS_STRATEGIES.get(bcs_class, BCS_STRATEGIES["I"])

    # Dosage form profile
    df_profile = DOSAGE_FORMS.get(dosage_form, DOSAGE_FORMS["immediate_release_tablet"])

    # Alternative dosage forms based on BCS
    alt_forms = []
    for form_key, form_prof in DOSAGE_FORMS.items():
        if bcs_class in form_prof.bcs_suitability and form_key != dosage_form:
            alt_forms.append(form_prof.form)

    # API-Excipient Compatibility Check
    compat_results = []
    for rule in INCOMPATIBILITY_RULES:
        for ing in blend:
            ing_l = ing.lower()
            # Check if ingredient matches trigger
            trigger_match = any(t.lower() in ing_l for t in rule.get("trigger_ingredients", []))
            func_match = False
            if ing in idx.index and "Function" in idx.columns:
                func = str(idx.loc[ing, "Function"])
                func_match = any(f.lower() in func.lower() for f in rule.get("trigger_functions", []))

            if trigger_match or func_match:
                compat_results.append(CompatibilityResult(
                    ingredient=ing,
                    compatibility_class="Conditional",
                    interaction_type=rule["interaction"],
                    mechanism=rule["mechanism"],
                    severity=rule["severity"],
                    mitigation=rule["mitigation"],
                    literature_ref=rule["ref"],
                ))

    # Compatibility summary
    severe = [c for c in compat_results if c.severity == "Severe"]
    moderate = [c for c in compat_results if c.severity == "Moderate"]
    if severe:
        compat_summary = f"⚠️ {len(severe)} severe compatibility concern(s) — review required"
    elif moderate:
        compat_summary = f"ℹ️ {len(moderate)} moderate compatibility note(s) — monitor"
    else:
        compat_summary = "✅ No significant compatibility concerns identified"

    # ICH Stability Zone
    primary_zone = "IVA" if any(m in target_markets for m in ["India", "Southeast Asia",
                                                                "Sub-Saharan Africa"]) else \
                   "II"  if any(m in target_markets for m in ["Southern Europe", "Japan"]) else "I"
    stability_profile = ICH_STABILITY_ZONES.get(primary_zone)

    stability_concerns = []
    packaging_rec = "HDPE bottle with desiccant (standard for oral solid dosage forms)"

    # Check blend for stability risk factors
    for ing in blend:
        ing_l = ing.lower()
        if "lactose" in ing_l:
            stability_concerns.append(f"{ing}: Reducing sugar — monitor for Maillard reaction with API")
        if "starch" in ing_l:
            stability_concerns.append(f"{ing}: Starch — moisture sorption may affect stability; "
                                      f"store below 25°C/60%RH")
        if "magnesium stearate" in ing_l:
            stability_concerns.append(f"{ing}: Metal stearate — may catalyze oxidation of "
                                      f"susceptible APIs at elevated temperature")
        if "peroxide" in ing_l or "percarbonate" in ing_l:
            stability_concerns.append(f"{ing}: Oxidizing agent — incompatible with "
                                      f"oxidation-sensitive APIs")

    if primary_zone in ["IVA", "IVB"]:
        packaging_rec = "Aluminium-Aluminium blister or HDPE with molecular sieve desiccant"
    elif primary_zone == "III":
        packaging_rec = "Hermetically sealed HDPE or glass with desiccant + antioxidant headspace"

    # Manufacturing route
    filler_pct = sum(pct for n, pct in blend.items() if n in idx.index and
                    any(x in str(idx.loc[n, "Function"] if "Function" in idx.columns else "")
                        for x in ["Filler", "Binder/Filler"]))
    binder_pct = sum(pct for n, pct in blend.items() if n in idx.index and
                    "Binder" in str(idx.loc[n, "Function"] if "Function" in idx.columns else ""))

    if filler_pct >= 50:
        mfg_route = "Direct Compression (DC)"
        mfg_rationale = (f"High filler content ({filler_pct:.1f}%) and MCC detected. "
                        f"DC is preferred — lower manufacturing cost, fewer steps, "
                        f"no thermal/moisture exposure.")
    elif binder_pct >= 5:
        mfg_route = "Wet Granulation"
        mfg_rationale = (f"Binder content ({binder_pct:.1f}%) suitable for wet granulation. "
                        f"Improves compressibility and content uniformity. "
                        f"Not recommended for moisture/heat-sensitive APIs.")
    else:
        mfg_route = "Dry Granulation / Roller Compaction"
        mfg_rationale = "Low filler content detected. Dry granulation avoids moisture exposure. "

    # Regulatory pathway
    if is_generic:
        pathway_key = "ANDA_505j"
    elif bcs_class == "I":
        pathway_key = "NDA_505b2"
    else:
        pathway_key = "NDA_505b1"
    pathway_profile = REGULATORY_PATHWAYS.get(pathway_key)

    # Development risks
    risks = []
    if bcs_class == "IV":
        risks.append("BCS Class IV — highest formulation challenge, consider alternative API form")
    if bcs_class == "II":
        risks.append("BCS Class II — solubility enhancement required, in vitro/in vivo correlation critical")
    if severe:
        risks.append(f"API-Excipient incompatibility: {', '.join(c.ingredient for c in severe)}")
    if is_pediatric:
        risks.append("Pediatric formulation: taste masking required, age-appropriate dosing critical")
    if primary_zone in ["IVA", "IVB"]:
        risks.append("Tropical stability zone: enhanced packaging and stability testing required")

    # Recommendations
    recommendations = bcs_profile.enabling_excipients[:3] + [
        f"Follow {mfg_route} manufacturing process",
        f"Conduct ICH Q1A stability at {stability_profile.long_term if stability_profile else '25°C/60%RH'}",
        f"Package in: {packaging_rec[:60]}",
    ]
    if compat_results:
        recommendations.append(f"Address {len(compat_results)} excipient compatibility note(s)")

    # Excipient recommendations based on dosage form
    exc_recs = []
    if "tablet" in dosage_form:
        exc_recs = [
            EXCIPIENT_RECOMMENDATIONS["filler_general"],
            EXCIPIENT_RECOMMENDATIONS["disintegrant"],
            EXCIPIENT_RECOMMENDATIONS["lubricant_standard"],
            EXCIPIENT_RECOMMENDATIONS["binder_wet_granulation"],
            EXCIPIENT_RECOMMENDATIONS["film_coat_ir"],
        ]
    elif "capsule" in dosage_form:
        exc_recs = [
            EXCIPIENT_RECOMMENDATIONS["filler_moisture_sensitive"],
            EXCIPIENT_RECOMMENDATIONS["lubricant_standard"],
        ]
    elif "solution" in dosage_form or "liquid" in dosage_form:
        exc_recs = [EXCIPIENT_RECOMMENDATIONS["preservative_liquid"]]
    elif "enteric" in dosage_form:
        exc_recs = [EXCIPIENT_RECOMMENDATIONS["enteric_coat"]]

    # ICM Score (0-100) — formulation complexity metric
    icm = 50.0
    icm += {"I": 0, "II": 20, "III": 15, "IV": 30}.get(bcs_class, 0)
    icm -= len([c for c in compat_results if c.severity == "Severe"]) * 10
    icm -= len(risks) * 3
    if is_pediatric: icm += 10  # more complex
    icm = max(0.0, min(100.0, icm))

    return PharmaDeepDiveResult(
        bcs_class=bcs_class,
        bcs_profile=bcs_profile,
        recommended_dosage_form=df_profile.form,
        dosage_form_profile=df_profile,
        alternative_forms=alt_forms[:3],
        excipient_recommendations=exc_recs,
        compatibility_results=compat_results,
        compatibility_summary=compat_summary,
        stability_zone=primary_zone,
        stability_profile=stability_profile,
        stability_concerns=stability_concerns,
        packaging_recommendation=packaging_rec,
        manufacturing_route=mfg_route,
        manufacturing_rationale=mfg_rationale,
        suggested_pathway=pathway_profile.pathway if pathway_profile else "NDA 505(b)(2)",
        pathway_profile=pathway_profile,
        development_risks=risks,
        development_recommendations=recommendations,
        icm_score=icm,
    )
