"""
modules/interaction_checker.py
Ingredient Interaction Safety Checker for IntelliForm.

Inspired by OpenMix's community interaction knowledge base
(github.com/vijayvkrishnan/openmix), this module checks blends
for known dangerous or counterproductive ingredient combinations
across personal care, food, pharmaceutical, industrial, and agricultural
verticals.

Coverage:
  - Oxidation/reduction conflicts
  - pH-driven degradation incompatibilities
  - Chelation conflicts
  - Emulsification incompatibilities
  - Regulatory co-restriction conflicts
  - Fermentation / microbial contamination risks

Severity levels:
  CRITICAL  — safety hazard (exothermic reaction, toxic byproduct)
  HIGH      — formulation failure (precipitation, phase separation, loss of efficacy)
  MEDIUM    — efficacy concern (competing mechanisms, partial neutralization)
  LOW       — regulatory note (co-restriction, label requirement)
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple


@dataclass
class Interaction:
    """A known incompatibility between two ingredient classes or specific ingredients."""
    ingredient_a_pattern: str   # substring/keyword match
    ingredient_b_pattern: str
    severity: str               # CRITICAL / HIGH / MEDIUM / LOW
    mechanism: str              # short chemistry explanation
    recommendation: str         # what to do
    vertical: Optional[str] = None  # None = all verticals


@dataclass
class InteractionFlag:
    """A triggered interaction for a specific blend."""
    ingredient_a: str
    ingredient_b: str
    severity: str
    mechanism: str
    recommendation: str


@dataclass
class InteractionReport:
    """Full interaction assessment for a blend."""
    flags: List[InteractionFlag]
    critical_count: int
    high_count: int
    medium_count: int
    low_count: int
    safe: bool                  # True if no CRITICAL or HIGH flags
    summary: str


# ── Interaction Database ───────────────────────────────────────────────────────
# Each entry: (pattern_a, pattern_b, severity, mechanism, recommendation, vertical)
# Patterns are case-insensitive substring matches against ingredient names.
# OpenMix-inspired + personal care chemistry + industrial + food safety literature.

_INTERACTIONS: List[Interaction] = [

    # ── OXIDIZER / REDUCER CONFLICTS ─────────────────────────────────────────
    Interaction(
        "benzoyl peroxide", "vitamin c",
        "CRITICAL",
        "Benzoyl peroxide oxidises ascorbic acid; generates free radicals and degrades both actives. "
        "Reaction can produce irritating byproducts.",
        "Do not co-formulate. Apply BPO and Vitamin C at different times of day.",
    ),
    Interaction(
        "benzoyl peroxide", "retinol",
        "CRITICAL",
        "BPO oxidises retinol to retinoic acid derivatives; unpredictable potency increase and irritation risk.",
        "Avoid direct combination. Physical separation in formulation required.",
    ),
    Interaction(
        "benzoyl peroxide", "antioxidant",
        "HIGH",
        "BPO oxidises any antioxidant, neutralising efficacy and generating reactive species.",
        "Remove antioxidants from BPO-containing formulations or encapsulate separately.",
    ),
    Interaction(
        "sodium hypochlorite", "ammonia",
        "CRITICAL",
        "Generates chloramine gas — toxic irritant and potential asphyxiant. Never mix.",
        "Never combine in industrial cleaning. Use separately with full rinse between applications.",
        vertical="industrial",
    ),
    Interaction(
        "sodium hypochlorite", "acid",
        "CRITICAL",
        "Bleach + acid releases chlorine gas (Cl₂) — toxic and corrosive. Acute inhalation hazard.",
        "Store separately. Never mix bleach with acidic cleaners. Full chemical segregation required.",
        vertical="industrial",
    ),
    Interaction(
        "hydrogen peroxide", "iron",
        "CRITICAL",
        "Fenton reaction: H₂O₂ + Fe²⁺ generates hydroxyl radicals (•OH) — highly destructive oxidant.",
        "Remove iron-containing ingredients or catalysts. Chelate trace iron with EDTA if necessary.",
    ),
    Interaction(
        "hydrogen peroxide", "copper peptide",
        "CRITICAL",
        "Copper ions catalyse H₂O₂ decomposition via Fenton mechanism, degrading both actives.",
        "Never combine peroxide actives with copper peptide complexes.",
        vertical="personal_care",
    ),
    Interaction(
        "peroxyacetic acid", "reducing agent",
        "CRITICAL",
        "Peroxyacid + reductant is an exothermic redox reaction with risk of ignition.",
        "Physical separation required. Do not co-formulate or mix in tank.",
        vertical="agricultural",
    ),

    # ── pH / ACID-BASE CONFLICTS ──────────────────────────────────────────────
    Interaction(
        "retinol", "aha",
        "HIGH",
        "AHAs (glycolic, lactic acid) create pH ≤3.5 which accelerates retinol oxidative degradation "
        "and increases conversion to retinoic acid — unpredictable irritation.",
        "Use retinol at pH 5.0–6.5. Apply AHA and retinol in separate steps.",
        vertical="personal_care",
    ),
    Interaction(
        "retinol", "bha",
        "HIGH",
        "BHA (salicylic acid) creates low-pH environment accelerating retinol degradation.",
        "Separate into distinct formulations or apply at different times.",
        vertical="personal_care",
    ),
    Interaction(
        "niacinamide", "vitamin c",
        "MEDIUM",
        "Ascorbic acid + niacinamide at elevated temperature can form niacin + nicotinamide — "
        "may cause transient flushing, though stability concern is overstated in modern literature.",
        "Stable at pH 5.5–7 and room temperature. Avoid high-temperature processing.",
        vertical="personal_care",
    ),
    Interaction(
        "sodium bicarbonate", "citric acid",
        "MEDIUM",
        "Acid-base effervescence reaction — generates CO₂, volume expansion, pH shift. "
        "Loss of both actives if premixed wet.",
        "Keep dry until point of use. Classic effervescent system — must be formulated as separate phases.",
    ),
    Interaction(
        "sodium bicarbonate", "lactic acid",
        "MEDIUM",
        "Acid-base neutralisation releases CO₂. pH buffering disrupts intended activity.",
        "pH adjustment after combining required; consider alternative pH system.",
    ),
    Interaction(
        "potassium sorbate", "ascorbic acid",
        "MEDIUM",
        "Ascorbic acid at low pH can reduce sorbate efficacy; combined oxidation-reduction products "
        "potential dehydrosorbate impurity formation.",
        "Monitor pH — maintain 3.5–5.5. Preferred combination in validated systems.",
        vertical="food",
    ),
    Interaction(
        "sodium benzoate", "ascorbic acid",
        "HIGH",
        "Benzene formation: sodium benzoate + ascorbic acid in aqueous system can generate benzene "
        "(IARC Group 1 carcinogen) — documented in soft drinks above 5 ppb.",
        "Avoid co-formulation. If required, use antioxidants other than ascorbic acid. "
        "Ensure pH < 3.5 is not combined with both ingredients.",
        vertical="food",
    ),
    Interaction(
        "lecithin", "strong acid",
        "MEDIUM",
        "Phospholipid hydrolysis accelerates at pH < 3. Lecithin emulsifier degrades, emulsion breaks.",
        "Maintain pH ≥ 4.5 in lecithin-stabilised emulsions.",
        vertical="food",
    ),

    # ── CHELATION / METAL CONFLICTS ───────────────────────────────────────────
    Interaction(
        "edta", "zinc",
        "MEDIUM",
        "EDTA chelates zinc ions, potentially reducing efficacy of zinc-based actives "
        "(e.g., zinc pyrithione in anti-dandruff, zinc oxide in sunscreen).",
        "Reduce EDTA concentration or use alternative chelant (GLDA) in zinc-active formulations.",
        vertical="personal_care",
    ),
    Interaction(
        "edta", "calcium",
        "LOW",
        "EDTA chelates calcium — beneficial for hard water but may strip mineralised substrates. "
        "Consider in oral care (dicalcium phosphate) applications.",
        "Evaluate EDTA level vs. calcium source concentration for intended application.",
    ),

    # ── EMULSIFICATION INCOMPATIBILITIES ─────────────────────────────────────
    Interaction(
        "cationic", "anionic",
        "HIGH",
        "Cationic and anionic surfactants form insoluble precipitates/complexes, "
        "phase separation, and loss of surface activity.",
        "Use non-ionic or zwitterionic bridging surfactants. Confirm HLB compatibility. "
        "Electrostatic incompatibility — do not co-formulate without compatibility testing.",
    ),
    Interaction(
        "quaternary ammonium", "anionic surfactant",
        "HIGH",
        "Quaternary ammonium compounds (quats) complex with anionic surfactants — "
        "precipitation, significant loss of antimicrobial and cleaning efficacy.",
        "Replace anionic surfactants with non-ionic alternatives when using quat actives.",
    ),

    # ── PRESERVATIVE CONFLICTS ────────────────────────────────────────────────
    Interaction(
        "parabens", "protein",
        "MEDIUM",
        "Parabens can interact with protein-based ingredients, reducing preservative efficacy "
        "through binding. More pronounced with high molecular weight proteins.",
        "Increase paraben concentration or switch to broad-spectrum alternatives. "
        "Perform challenge testing per ISO 11930.",
        vertical="personal_care",
    ),
    Interaction(
        "phenoxyethanol", "high peg",
        "MEDIUM",
        "High PEG/polyol concentrations solubilise phenoxyethanol away from the aqueous phase, "
        "reducing antimicrobial efficacy (Minimum Inhibitory Concentration increases).",
        "Use phenoxyethanol at higher %. Perform preservative efficacy testing when PEG > 5%.",
    ),

    # ── OXIDATION / STABILITY CONFLICTS ──────────────────────────────────────
    Interaction(
        "sulfur dust", "organic solvent",
        "HIGH",
        "Elemental sulfur is incompatible with many organic solvents at elevated temperatures — "
        "exothermic reactions possible. Sulfur + strong oxidizers risk of ignition.",
        "Use aqueous slurry or encapsulated sulfur. Avoid mixing with volatile organic carriers.",
        vertical="agricultural",
    ),
    Interaction(
        "sodium percarbonate", "enzyme",
        "HIGH",
        "Percarbonate releases H₂O₂ which denatures proteases, amylases, and lipases. "
        "Enzyme activity rapidly destroyed in presence of active oxygen bleach.",
        "Add enzymes after percarbonate has decomposed, or use encapsulated enzymes.",
        vertical="fabric_laundry",
    ),
    Interaction(
        "sodium percarbonate", "reducing agent",
        "HIGH",
        "Oxidation-reduction reaction neutralises both actives; exothermic risk at high concentrations.",
        "Do not combine with reducing agents in concentrated form.",
    ),

    # ── AGRICULTURAL SPECIFIC ─────────────────────────────────────────────────
    Interaction(
        "calcium hydroxide", "sulfur",
        "MEDIUM",
        "Lime-sulfur (calcium polysulfide) reaction — forms phytotoxic polysulfides at high temperature. "
        "Known crop burn risk if applied together.",
        "Do not tank-mix lime and sulfur. Apply separately with 48-hour interval.",
        vertical="agricultural",
    ),
    Interaction(
        "copper sulfate", "phosphate",
        "MEDIUM",
        "Copper-phosphate precipitation reduces bioavailability of both copper and phosphate. "
        "Reduced efficacy of copper-based fungicides and fertiliser phosphorus.",
        "Separate application timing or use chelated copper formulations.",
        vertical="agricultural",
    ),

    # ── INDUSTRIAL SPECIFIC ───────────────────────────────────────────────────
    Interaction(
        "strong alkali", "aluminum",
        "HIGH",
        "NaOH/KOH react with aluminium to generate hydrogen gas (H₂) — explosion risk in confined spaces.",
        "Avoid contact of strong caustic with aluminium equipment or surfaces. "
        "Use plastic or stainless steel equipment.",
        vertical="industrial",
    ),
    Interaction(
        "oxidizing acid", "flammable solvent",
        "CRITICAL",
        "Nitric/chromic acids with organic solvents — highly exothermic, fire and explosion risk.",
        "Strict physical separation. Never mix oxidising acids with organic solvents.",
        vertical="industrial",
    ),
]


def _matches(ingredient: str, pattern: str) -> bool:
    """Case-insensitive substring match of pattern against ingredient name."""
    return pattern.lower() in ingredient.lower()


def check_interactions(
    blend: dict,
    vertical: Optional[str] = None,
) -> InteractionReport:
    """
    Check all ingredient pairs in a blend for known incompatibilities.

    Args:
        blend: {ingredient_name: weight_percent}
        vertical: optional vertical key to filter vertical-specific rules

    Returns:
        InteractionReport with all triggered flags.
    """
    ingredients = list(blend.keys())
    flags: List[InteractionFlag] = []

    for rule in _INTERACTIONS:
        # Skip vertical-specific rules that don't match
        if rule.vertical and vertical and rule.vertical != vertical:
            continue

        # Check all ingredient pairs against this rule
        matched_a: List[str] = [i for i in ingredients if _matches(i, rule.ingredient_a_pattern)]
        matched_b: List[str] = [i for i in ingredients if _matches(i, rule.ingredient_b_pattern)]

        if matched_a and matched_b:
            for a in matched_a:
                for b in matched_b:
                    if a != b:
                        flags.append(InteractionFlag(
                            ingredient_a=a,
                            ingredient_b=b,
                            severity=rule.severity,
                            mechanism=rule.mechanism,
                            recommendation=rule.recommendation,
                        ))

    # Deduplicate
    seen: Set[Tuple[str, str, str]] = set()
    unique_flags: List[InteractionFlag] = []
    for f in flags:
        key = (min(f.ingredient_a, f.ingredient_b), max(f.ingredient_a, f.ingredient_b), f.severity)
        if key not in seen:
            seen.add(key)
            unique_flags.append(f)

    critical = sum(1 for f in unique_flags if f.severity == "CRITICAL")
    high     = sum(1 for f in unique_flags if f.severity == "HIGH")
    medium   = sum(1 for f in unique_flags if f.severity == "MEDIUM")
    low      = sum(1 for f in unique_flags if f.severity == "LOW")

    safe = critical == 0 and high == 0

    if critical > 0:
        summary = f"⛔ {critical} CRITICAL interaction(s) detected — blend is NOT SAFE for use. Review immediately."
    elif high > 0:
        summary = f"⚠️  {high} HIGH-severity interaction(s) — formulation failure or significant efficacy loss expected."
    elif medium > 0:
        summary = f"⚡ {medium} MEDIUM concern(s) — efficacy or stability may be impacted."
    elif low > 0:
        summary = f"ℹ️  {low} LOW-severity note(s) — regulatory or label awareness."
    else:
        summary = "✅ No known ingredient interactions detected."

    return InteractionReport(
        flags=unique_flags,
        critical_count=critical,
        high_count=high,
        medium_count=medium,
        low_count=low,
        safe=safe,
        summary=summary,
    )
