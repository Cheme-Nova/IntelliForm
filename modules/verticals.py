"""
modules/verticals.py
Vertical-specific constraint profiles for IntelliForm v1.0.

Each vertical has:
  - display name and icon
  - default optimizer constraints (cost, min_bio, min_perf)
  - regulatory framework tags
  - required function types
  - excluded function types
  - specific constraint flags
"""
from dataclasses import dataclass, field
from typing import List, Dict, Optional


@dataclass
class VerticalProfile:
    key: str
    name: str
    icon: str
    description: str

    # Default optimizer constraints
    default_max_cost: float
    default_min_bio: float
    default_min_perf: float

    # Regulatory frameworks to check
    regulatory_frameworks: List[str]

    # Required ingredient functions (at least one must be present)
    required_functions: List[str]

    # Functions to exclude from this vertical
    excluded_functions: List[str]

    # Specific constraint flags
    constraints: Dict[str, str] = field(default_factory=dict)

    # Example prompts for the sidebar
    example_prompts: List[str] = field(default_factory=list)


VERTICAL_PROFILES = {

    "personal_care": VerticalProfile(
        key="personal_care",
        name="Personal Care & Cosmetics",
        icon="💄",
        description="Shampoos, conditioners, cleansers, lotions, serums",
        default_max_cost=6.0,
        default_min_bio=95.0,
        default_min_perf=82.0,
        regulatory_frameworks=["REACH", "EU_ECOLABEL", "COSMOS", "EPA_SAFER_CHOICE"],
        required_functions=["Surfactant", "Emulsifier", "Humectant", "Emollient",
                           "Thickener", "Preservative", "Conditioning"],
        excluded_functions=["Pharma", "Agri", "Paint", "Fabric Stain",
                           "Biopesticide", "Bioinsecticide"],
        constraints={
            "skin_safe": "All ingredients must be skin-safe",
            "dermatologically_tested": "Mild formulation priority",
        },
        example_prompts=[
            "Mild foaming cleanser for sensitive skin, 98% bio-based, under $5/kg",
            "Moisturizing body lotion with hyaluronic acid, COSMOS certified",
            "Anti-dandruff shampoo with natural actives, sulfate-free",
        ]
    ),

    "industrial": VerticalProfile(
        key="industrial",
        name="Industrial Cleaning",
        icon="🏭",
        description="Degreasers, metal cleaners, hard surface cleaners, disinfectants",
        default_max_cost=4.0,
        default_min_bio=80.0,
        default_min_perf=85.0,
        regulatory_frameworks=["REACH", "EPA_SAFER_CHOICE", "NSF_A1"],
        required_functions=["Industrial Solvent", "Industrial Cleaner", "Surfactant",
                           "Chelating Agent", "pH Adjuster", "Degreaser",
                           "Bio-Solvent", "Hydrotrope"],
        excluded_functions=["Pharma", "Food", "Agri Biopesticide", "Paint Binder",
                           "Skin Active", "Skin Conditioning", "Emollient"],
        constraints={
            "drain_safe": "Municipal drain discharge compliance",
            "low_voc": "Low VOC formulation preferred",
            "worker_safe": "Worker safety — low inhalation risk",
        },
        example_prompts=[
            "Heavy duty metal degreaser, EPA Safer Choice, drain-safe, under $3.50/kg",
            "Food-safe kitchen degreaser, NSF A1 approved, 90% bio-based",
            "Concrete cleaner for industrial floors, phosphate-free, low foam",
        ]
    ),

    "agricultural": VerticalProfile(
        key="agricultural",
        name="Agricultural",
        icon="🌾",
        description="Crop protection, adjuvants, biostimulants, fertilizers",
        default_max_cost=5.0,
        default_min_bio=85.0,
        default_min_perf=80.0,
        regulatory_frameworks=["EPA_FIFRA", "EU_REGULATION_1107", "OMRI"],
        required_functions=["Agri", "Agricultural", "Biopesticide", "Biostimulant",
                           "Adjuvant", "Bioinsecticide", "Biocontrol"],
        excluded_functions=["Pharma Active", "Food Sweetener", "Paint Binder",
                           "Skin Active", "Hair Conditioning", "Fabric Softener"],
        constraints={
            "rei_compliant": "Restricted Entry Interval compliance",
            "phi_compliant": "Pre-Harvest Interval compliance",
            "bee_safe": "Low toxicity to pollinators",
            "soil_safe": "Minimal soil persistence",
        },
        example_prompts=[
            "OMRI-listed bioinsecticide adjuvant for organic crops, under $8/kg",
            "Biostimulant blend for root growth, 100% bio-based, soil-safe",
            "Spreader-sticker for herbicide tank mix, drift retardant",
        ]
    ),

    "pharmaceutical": VerticalProfile(
        key="pharmaceutical",
        name="Pharmaceutical",
        icon="💊",
        description="Tablet excipients, capsule fillers, topical bases, oral solutions",
        default_max_cost=15.0,
        default_min_bio=80.0,
        default_min_perf=78.0,
        regulatory_frameworks=["ICH_Q8", "USP_NF", "Ph_Eur", "FDA_21CFR"],
        required_functions=["Pharma", "Excipient", "Binder", "Filler", "Disintegrant",
                           "Lubricant", "Film Coat", "Solubilizer", "Buffer"],
        excluded_functions=["Agri Biopesticide", "Industrial Degreaser", "Paint Binder",
                           "Fabric Softener", "Optical Brightener"],
        constraints={
            "gmp_compliant": "GMP grade materials only",
            "compendial": "USP/NF or Ph.Eur compendial preferred",
            "api_compatible": "Excipient-API compatibility required",
            "endotoxin": "Low endotoxin for parenteral routes",
        },
        example_prompts=[
            "Direct compression tablet blend, 500mg tablet, MCC + disintegrant + lubricant",
            "Oral liquid suspension, pH 4-6, preservative system, USP grade",
            "Enteric-coated pellet formulation, acid-resistant, sustained release",
        ]
    ),

    "food": VerticalProfile(
        key="food",
        name="Food & Beverage",
        icon="🍃",
        description="Food additives, functional ingredients, flavor systems, preservatives",
        default_max_cost=5.0,
        default_min_bio=90.0,
        default_min_perf=75.0,
        regulatory_frameworks=["FDA_GRAS", "EU_FOOD_ADDITIVE", "EFSA", "CODEX"],
        required_functions=["Food", "GRAS", "Emulsifier", "Thickener", "Stabilizer",
                           "Acidulant", "Preservative", "Sweetener", "Flavor",
                           "Color", "Protein", "Fiber"],
        excluded_functions=["Pharma Active", "Industrial Degreaser", "Paint Binder",
                           "Agri Biopesticide", "Fabric Softener",
                           "Industrial Corrosion", "Industrial Biocide"],
        constraints={
            "gras_status": "GRAS or approved food additive status",
            "allergen_free": "Major allergen declaration required",
            "clean_label": "Clean label ingredients preferred",
            "halal_kosher": "Halal/Kosher certifiable preferred",
        },
        example_prompts=[
            "Clean-label emulsifier blend for plant-based milk, GRAS, under $4/kg",
            "Natural preservative system for bakery, no artificial additives",
            "Fiber-enriched thickener for functional beverage, neutral flavor",
        ]
    ),

    "fabric_laundry": VerticalProfile(
        key="fabric_laundry",
        name="Fabric & Laundry",
        icon="👕",
        description="Detergents, fabric softeners, stain removers, laundry additives",
        default_max_cost=3.5,
        default_min_bio=85.0,
        default_min_perf=82.0,
        regulatory_frameworks=["REACH", "EU_DETERGENT_REGULATION", "EPA_DESIGN_FOR_ENVIRONMENT"],
        required_functions=["Fabric", "Surfactant", "Builder", "Enzyme", "Softener",
                           "Bleach", "Optical Brightener", "Anti-Redeposition",
                           "Soil Release", "Fragrance"],
        excluded_functions=["Pharma Active", "Food Sweetener", "Paint Binder",
                           "Skin Active", "Agri Biopesticide"],
        constraints={
            "readily_biodegradable": ">60% biodegradation in 28 days",
            "phosphate_free": "No phosphate builders",
            "low_temp_active": "Active at 30-40°C wash temperatures",
            "color_safe": "Color fabric safe",
        },
        example_prompts=[
            "Concentrated liquid detergent, phosphate-free, cold-water active, 85% bio-based",
            "Plant-based fabric softener with esterquat, biodegradable, floral fragrance",
            "Enzyme-boosted stain remover, protease + lipase + amylase, under $3/kg",
        ]
    ),

    "paint_coatings": VerticalProfile(
        key="paint_coatings",
        name="Paint & Coatings",
        icon="🎨",
        description="Waterborne paints, varnishes, industrial coatings, wood finishes",
        default_max_cost=6.0,
        default_min_bio=60.0,
        default_min_perf=80.0,
        regulatory_frameworks=["REACH", "EU_DECOPAINT_DIRECTIVE", "VOC_REGULATION", "LEED"],
        required_functions=["Paint", "Coatings", "Binder", "Pigment", "Extender",
                           "Rheology Modifier", "Coalescent", "Dispersant",
                           "Defoamer", "UV Stabilizer", "Crosslinker"],
        excluded_functions=["Pharma Active", "Food Sweetener", "Agri Biopesticide",
                           "Fabric Softener", "Skin Active"],
        constraints={
            "low_voc": "VOC <50g/L (decorative) or <250g/L (industrial)",
            "biocide_free": "Dry film biocide restrictions",
            "freeze_thaw": "Freeze-thaw stability",
            "scrub_resistance": "Wet scrub resistance >200 cycles",
        },
        example_prompts=[
            "Zero-VOC waterborne wall paint binder, high scrub resistance, LEED compliant",
            "Bio-based alkyd varnish for wood, outdoor durable, low yellowing",
            "Industrial protective coating, corrosion resistant, bio-based content 60%",
        ]
    ),

}

# All vertical keys and display names for selector
VERTICAL_OPTIONS = {
    "all": "🌐 All Verticals",
    "personal_care": "💄 Personal Care & Cosmetics",
    "industrial": "🏭 Industrial Cleaning",
    "agricultural": "🌾 Agricultural",
    "pharmaceutical": "💊 Pharmaceutical",
    "food": "🍃 Food & Beverage",
    "fabric_laundry": "👕 Fabric & Laundry",
    "paint_coatings": "🎨 Paint & Coatings",
}


def get_profile(vertical_key: str) -> Optional[VerticalProfile]:
    return VERTICAL_PROFILES.get(vertical_key)


def filter_db_by_vertical(db, vertical_key: str):
    """Filter ingredient database to only show relevant vertical."""
    import pandas as pd
    if vertical_key == "all" or vertical_key not in VERTICAL_PROFILES:
        return db
    if "Vertical" not in db.columns:
        return db

    # Primary: ingredients tagged for this vertical
    filtered = db[db["Vertical"] == vertical_key].copy()

    # Secondary: only pull truly universal ingredients — water, basic salts, citric acid
    # Do NOT include industrial solvents, pharma vehicles, etc. as "common"
    UNIVERSAL = [
        "Water", "Citric Acid", "Lactic Acid", "Sodium Chloride",
        "Glycerol", "Sodium Bicarbonate", "Sodium Citrate",
        "Potassium Citrate", "Malic Acid", "Tartaric Acid",
    ]
    universal_mask = db["Ingredient"].apply(
        lambda x: any(u.lower() in x.lower() for u in UNIVERSAL)
    )
    universal = db[universal_mask]

    combined = pd.concat([filtered, universal]).drop_duplicates(subset=["Ingredient"])
    return combined.reset_index(drop=True)


def get_vertical_constraints(vertical_key: str, parsed_cost, parsed_bio, parsed_perf):
    """
    Override parsed constraints with vertical minimums where vertical is stricter.
    Returns (max_cost, min_bio, min_perf)
    """
    profile = get_profile(vertical_key)
    if not profile:
        return parsed_cost, parsed_bio, parsed_perf
    # Use vertical defaults if parsed values are more permissive
    max_cost = min(parsed_cost, profile.default_max_cost) if parsed_cost < profile.default_max_cost * 1.5 else parsed_cost
    min_bio  = max(parsed_bio,  profile.default_min_bio)
    min_perf = max(parsed_perf, profile.default_min_perf)
    return max_cost, min_bio, min_perf
