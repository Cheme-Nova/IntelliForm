"""
modules/sdl_integration.py
Self-Driving Lab (SDL) Integration — IntelliForm v1.6

Bridges IntelliForm formulation outputs to autonomous/automated lab hardware
and closes the loop with Closed-Loop Reformulation Intelligence
(modules/reformulation_intelligence.py).

What it does:
  1. generate_sdl_protocol() — turns a {ingredient: %} blend into a
     machine-readable liquid-handling protocol (generic JSON + an
     Opentrons OT-2 Python protocol stub) for a chosen batch scale.
  2. ingest_sdl_results() — takes measured properties back from the SDL
     (or a manual pilot batch QC sheet), compares against target specs,
     auto-classifies any failure using the same taxonomy as
     reformulation_intelligence, and — if a failure is detected — runs
     the reformulation engine to propose the next blend to test.
  3. run_closed_loop() — drives propose → run → measure → reformulate →
     repeat for up to N iterations, given a `measurement_fn` callable
     (the real SDL driver, or a simulator for demos).

This module does not talk to lab hardware directly — it produces and
consumes structured data so it can sit behind any SDL orchestration layer
(Opentrons HTTP API, Chemspeed, a manual lab notebook, etc.).

Reference:
  - Self-driving labs: Abolhasani & Kumacheva, Nature Synthesis, 2023
  - Opentrons Python Protocol API v2: https://docs.opentrons.com
"""
from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime
from typing import Callable, Dict, List, Optional, Tuple
import uuid

import pandas as pd

from modules.reformulation_intelligence import (
    run_reformulation_intelligence, FAILURE_TYPES, ReformulationReport,
)
from modules.analytics import track


SUPPORTED_PLATFORMS = ["opentrons_ot2", "generic_liquid_handler", "manual"]

# Coarse density hints (g/mL) so mass % can be converted to dispense volumes.
# Defaults to 1.0 (water-like) when no keyword matches.
DENSITY_HINTS: List[Tuple[str, float]] = [
    ("oil", 0.92), ("butter", 0.93), ("wax", 0.95), ("ester", 0.88),
    ("silicone", 0.96), ("glycerin", 1.26), ("glycerol", 1.26),
    ("ethanol", 0.79), ("alcohol", 0.85),
]

# Maps measured-property keys → (failure_type if too high, failure_type if too low)
SPEC_FAILURE_MAP: Dict[str, Tuple[str, str]] = {
    "ph":          ("ph_too_high", "ph_too_low"),
    "viscosity_cP": ("viscosity_too_high", "viscosity_too_low"),
}


def _density_for(ingredient: str) -> float:
    ing_l = ingredient.lower()
    for kw, density in DENSITY_HINTS:
        if kw in ing_l:
            return density
    return 1.0


# ── Result schemas ───────────────────────────────────────────────────────────

@dataclass
class LiquidTransferStep:
    step: int
    ingredient: str
    target_mass_g: float
    target_volume_uL: float
    density_g_mL: float
    source_labware: str
    dest_labware: str
    notes: str


@dataclass
class SDLProtocol:
    protocol_id: str
    platform: str
    batch_scale_g: float
    blend: Dict[str, float]
    transfers: List[LiquidTransferStep]
    measurement_requests: List[str]
    estimated_runtime_min: float
    opentrons_snippet: Optional[str]
    created_at: str


@dataclass
class SDLLoopResult:
    iteration: int
    measured: Dict[str, float]
    failure_detected: Optional[str]
    failure_label: Optional[str]
    reformulation: Optional[ReformulationReport]
    next_blend: Dict[str, float]
    status: str           # "pass" | "reformulate" | "needs_review"
    notes: str


# ── Protocol generation ──────────────────────────────────────────────────────

def generate_sdl_protocol(
    blend: Dict[str, float],
    db: pd.DataFrame,
    batch_scale_g: float = 10.0,
    platform: str = "opentrons_ot2",
    target_specs: Optional[Dict[str, float]] = None,
) -> SDLProtocol:
    """
    Generate an SDL-ready experiment protocol for a formulation blend.

    Args:
        blend: {ingredient: percentage} — must sum to ~100
        db: ingredient database (currently unused, reserved for stock-concentration lookups)
        batch_scale_g: total mass of the experiment in grams (e.g. 10 g well-plate scale)
        platform: one of SUPPORTED_PLATFORMS
        target_specs: optional dict of target specs (e.g. {"target_ph_min": 5, ...})
            used to populate measurement_requests

    Returns:
        SDLProtocol with per-ingredient liquid transfer steps and, for
        opentrons_ot2, a runnable Python protocol stub.
    """
    if platform not in SUPPORTED_PLATFORMS:
        platform = "generic_liquid_handler"

    total_pct = sum(blend.values()) or 100.0
    transfers: List[LiquidTransferStep] = []
    for i, (ing, pct) in enumerate(blend.items(), start=1):
        mass_g = batch_scale_g * (pct / total_pct)
        density = _density_for(ing)
        volume_uL = (mass_g / density) * 1000.0
        transfers.append(LiquidTransferStep(
            step=i, ingredient=ing,
            target_mass_g=round(mass_g, 4),
            target_volume_uL=round(volume_uL, 1),
            density_g_mL=density,
            source_labware="stock_reservoir",
            dest_labware="reaction_vial_1",
            notes="Volume estimated from mass % using density heuristic — "
                  "calibrate against actual stock density before SDL run.",
        ))

    # Measurement requests derived from target specs (fallback to a standard panel)
    measurement_requests: List[str] = []
    if target_specs:
        if any(k.startswith("target_ph") for k in target_specs):
            measurement_requests.append("pH (probe, 25°C)")
        if any("viscosity" in k for k in target_specs):
            measurement_requests.append("Viscosity (cP, spindle/shear rate per SOP)")
        if any("performance" in k for k in target_specs):
            measurement_requests.append(f"Performance: {target_specs.get('performance_metric', 'primary KPI')}")
        if any("colour" in k or "color" in k for k in target_specs):
            measurement_requests.append("Appearance / colour (visual or spectrophotometer)")
    if not measurement_requests:
        measurement_requests = ["pH (probe, 25°C)", "Viscosity (cP)", "Appearance (visual)"]

    # Runtime estimate: ~1.5 min per transfer + 5 min mix/measure overhead
    runtime_min = round(len(transfers) * 1.5 + 5.0, 1)

    protocol_id = f"SDL-{datetime.now().strftime('%Y%m%d%H%M%S')}-{uuid.uuid4().hex[:6]}"

    opentrons_snippet = None
    if platform == "opentrons_ot2":
        opentrons_snippet = _build_opentrons_snippet(protocol_id, transfers, batch_scale_g)

    track("sdl_protocol_generated", {
        "platform": platform, "batch_scale_g": batch_scale_g,
        "n_transfers": len(transfers),
    })

    return SDLProtocol(
        protocol_id=protocol_id,
        platform=platform,
        batch_scale_g=batch_scale_g,
        blend=blend,
        transfers=transfers,
        measurement_requests=measurement_requests,
        estimated_runtime_min=runtime_min,
        opentrons_snippet=opentrons_snippet,
        created_at=datetime.now().isoformat(),
    )


def _build_opentrons_snippet(protocol_id: str, transfers: List[LiquidTransferStep],
                              batch_scale_g: float) -> str:
    """Build a minimal Opentrons Python Protocol API v2 stub for review/customization."""
    lines = [
        '"""',
        f"IntelliForm SDL Protocol — {protocol_id}",
        f"Batch scale: {batch_scale_g} g — generated by modules/sdl_integration.py",
        "Review labware/slot assignments before running on hardware.",
        '"""',
        "from opentrons import protocol_api",
        "",
        "metadata = {",
        f'    "protocolName": "{protocol_id}",',
        '    "author": "IntelliForm SDL Integration",',
        '    "apiLevel": "2.18",',
        "}",
        "",
        "def run(protocol: protocol_api.ProtocolContext):",
        '    tiprack = protocol.load_labware("opentrons_96_tiprack_300ul", 1)',
        '    reservoir = protocol.load_labware("nest_12_reservoir_15ml", 2)',
        '    dest = protocol.load_labware("nest_96_wellplate_200ul_flat", 3)',
        '    pipette = protocol.load_instrument("p300_single_gen2", "right", tip_racks=[tiprack])',
        "",
        "    pipette.pick_up_tip()",
    ]
    for t in transfers:
        lines.append(
            f'    # {t.ingredient}: {t.target_mass_g} g (~{t.target_volume_uL} uL, '
            f'rho={t.density_g_mL} g/mL)'
        )
        lines.append(
            f'    pipette.transfer({t.target_volume_uL}, reservoir["A1"], dest["A1"], new_tip="never")'
        )
    lines.append("    pipette.drop_tip()")
    return "\n".join(lines)


# ── Result ingestion / closed loop ───────────────────────────────────────────

def _detect_failure(measured: Dict[str, float],
                     target_specs: Dict[str, float]) -> Tuple[Optional[str], dict]:
    """
    Compare measured properties against target specs and classify any failure
    using the reformulation_intelligence FAILURE_TYPES taxonomy.

    Returns (failure_type | None, test_data) where test_data is ready to pass
    to run_reformulation_intelligence().
    """
    # pH check
    if "measured_ph" in measured or "ph" in measured:
        ph = measured.get("measured_ph", measured.get("ph"))
        ph_min = target_specs.get("target_ph_min")
        ph_max = target_specs.get("target_ph_max")
        if ph is not None and ph_min is not None and ph_max is not None:
            if ph > ph_max:
                return "ph_too_high", {
                    "measured_ph": ph, "target_ph_min": ph_min, "target_ph_max": ph_max,
                }
            if ph < ph_min:
                return "ph_too_low", {
                    "measured_ph": ph, "target_ph_min": ph_min, "target_ph_max": ph_max,
                }

    # Viscosity check
    visc = measured.get("measured_viscosity_cP", measured.get("viscosity_cP"))
    visc_max = target_specs.get("target_viscosity_max_cP")
    visc_min = target_specs.get("target_viscosity_min_cP")
    if visc is not None:
        if visc_max is not None and visc > visc_max:
            return "viscosity_too_high", {
                "measured_viscosity_cP": visc, "target_viscosity_max_cP": visc_max,
            }
        if visc_min is not None and visc < visc_min:
            return "viscosity_too_low", {
                "measured_viscosity_cP": visc, "target_viscosity_min_cP": visc_min,
            }

    # Phase separation flag (boolean from SDL imaging/turbidity check)
    if measured.get("phase_separation"):
        return "phase_separation", {
            "storage_days": measured.get("storage_days", 1),
            "temperature_C": measured.get("temperature_C", 45),
        }

    # Performance shortfall
    perf = measured.get("measured_performance")
    target_perf = target_specs.get("target_performance")
    if perf is not None and target_perf is not None and perf < target_perf:
        return "performance_shortfall", {
            "measured_performance": perf, "target_performance": target_perf,
            "performance_metric": target_specs.get("performance_metric", "primary KPI"),
        }

    return None, {}


def _apply_best_suggestion(blend: Dict[str, float],
                            report: ReformulationReport) -> Dict[str, float]:
    """Apply the top-ranked reformulation suggestion to produce the next blend to test."""
    next_blend = dict(blend)
    sug = report.best_suggestion
    if sug.action_type in ("increase", "decrease", "replace") and sug.ingredient in next_blend:
        next_blend[sug.ingredient] = round(sug.suggested_pct, 3)
    elif sug.action_type == "add" and sug.suggested_pct > 0:
        next_blend[sug.ingredient] = round(sug.suggested_pct, 3)
    elif sug.action_type == "remove":
        next_blend.pop(sug.ingredient, None)

    # Renormalize to 100%
    total = sum(next_blend.values())
    if total > 0:
        next_blend = {k: round(v * 100.0 / total, 3) for k, v in next_blend.items()}
    return next_blend


def ingest_sdl_results(
    blend: Dict[str, float],
    db: pd.DataFrame,
    measured: Dict[str, float],
    target_specs: Dict[str, float],
    vertical: str = "all",
    iteration: int = 1,
    batch_id: Optional[str] = None,
) -> SDLLoopResult:
    """
    Ingest a single SDL/pilot-batch measurement result, auto-detect any
    specification failure, and (if found) run Closed-Loop Reformulation
    Intelligence to propose the next blend to test.
    """
    failure_type, test_data = _detect_failure(measured, target_specs)

    if failure_type is None:
        track("sdl_result_pass", {"iteration": iteration, "vertical": vertical})
        return SDLLoopResult(
            iteration=iteration, measured=measured, failure_detected=None,
            failure_label=None, reformulation=None, next_blend=dict(blend),
            status="pass",
            notes="All measured properties within target specs — no reformulation needed.",
        )

    if failure_type not in FAILURE_TYPES:
        return SDLLoopResult(
            iteration=iteration, measured=measured, failure_detected=failure_type,
            failure_label=failure_type, reformulation=None, next_blend=dict(blend),
            status="needs_review",
            notes=f"Unrecognized failure signature '{failure_type}' — manual review required.",
        )

    report = run_reformulation_intelligence(
        blend=blend, db=db, failure_type=failure_type,
        test_data=test_data, batch_id=batch_id,
    )
    next_blend = _apply_best_suggestion(blend, report)

    track("sdl_result_reformulate", {
        "iteration": iteration, "failure_type": failure_type, "vertical": vertical,
    })

    return SDLLoopResult(
        iteration=iteration, measured=measured, failure_detected=failure_type,
        failure_label=FAILURE_TYPES[failure_type]["label"],
        reformulation=report, next_blend=next_blend,
        status="reformulate",
        notes=f"{FAILURE_TYPES[failure_type]['label']} detected — "
              f"applied top suggestion ({report.best_suggestion.action_type} "
              f"{report.best_suggestion.ingredient}) to produce next blend.",
    )


def run_closed_loop(
    blend: Dict[str, float],
    db: pd.DataFrame,
    target_specs: Dict[str, float],
    measurement_fn: Callable[[Dict[str, float], int], Dict[str, float]],
    vertical: str = "all",
    max_iterations: int = 5,
    batch_id_prefix: Optional[str] = None,
) -> List[SDLLoopResult]:
    """
    Drive the propose → run → measure → reformulate loop for up to
    `max_iterations` rounds.

    Args:
        blend: starting formulation {ingredient: %}
        db: ingredient database
        target_specs: target specification dict (see _detect_failure)
        measurement_fn: callable(blend, iteration) -> measured dict.
            In production this wraps the real SDL/instrument driver; for
            demos pass a simulator that perturbs results based on the blend.
        max_iterations: stop after this many rounds even if not passing
        batch_id_prefix: optional prefix for generated batch IDs

    Returns:
        List of SDLLoopResult, one per iteration (stops early on "pass").
    """
    history: List[SDLLoopResult] = []
    current_blend = dict(blend)

    for i in range(1, max_iterations + 1):
        measured = measurement_fn(current_blend, i)
        batch_id = f"{batch_id_prefix or 'SDL'}-{i:02d}"
        result = ingest_sdl_results(
            current_blend, db, measured, target_specs,
            vertical=vertical, iteration=i, batch_id=batch_id,
        )
        history.append(result)

        if result.status == "pass":
            break
        if result.status == "needs_review":
            break
        current_blend = result.next_blend

    track("sdl_closed_loop_complete", {
        "iterations_run": len(history),
        "final_status": history[-1].status if history else "none",
        "vertical": vertical,
    })
    return history
