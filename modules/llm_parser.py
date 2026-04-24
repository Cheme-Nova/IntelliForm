"""
modules/llm_parser.py
Natural language -> formulation constraints parser for IntelliForm.

This parser borrows the reliability lessons from FormulAI:
  - avoid brittle provider-side forced JSON modes
  - extract and repair JSON locally when models wrap or slightly damage it
  - normalize application types to IntelliForm's canonical vertical keys
  - validate and clamp numeric outputs before the optimizer sees them
"""
import json
import os
import re
import traceback
from dataclasses import dataclass
from typing import Optional

import requests

from modules.analytics import track
from modules.verticals import get_profile


@dataclass
class ParseResult:
    max_cost: float
    min_bio: float
    min_perf: float
    application_type: str
    reasoning: str
    parser_backend: str
    raw_input: str


VALID_APPLICATION_TYPES = {
    "personal_care",
    "industrial",
    "agricultural",
    "pharmaceutical",
    "food",
    "fabric_laundry",
    "paint_coatings",
    "unknown",
}


APPLICATION_TYPE_ALIASES = {
    "cosmetics": "personal_care",
    "cosmetic": "personal_care",
    "personal care": "personal_care",
    "home_care": "fabric_laundry",
    "home care": "fabric_laundry",
    "laundry": "fabric_laundry",
    "textile": "fabric_laundry",
    "textiles": "fabric_laundry",
    "pharma": "pharmaceutical",
    "food_safe": "food",
    "food safe": "food",
    "coatings": "paint_coatings",
    "paint": "paint_coatings",
}


_SYSTEM_PROMPT = """You are IntelliForm, a formulation-intelligence parser for specialty chemical R&D.
Extract optimization constraints from a customer brief.
Return exactly one JSON object and nothing else.
Do not use markdown, code fences, or explanatory text outside the JSON.

JSON schema:
{
  "max_cost": <float, USD/kg, default 5.0>,
  "min_bio": <float, bio-based percentage 0-100, default 95>,
  "min_perf": <float, performance score 0-100, default 82>,
  "application_type": <string, one of: personal_care | industrial | agricultural | pharmaceutical | food | fabric_laundry | paint_coatings | unknown>,
  "reasoning": <string, 1-2 short sentences>
}

Rules:
- Never set min_bio below 80 unless the brief is clearly pharmaceutical or explicitly low-bio.
- Never set min_perf below 65.
- Extract cost ceilings from phrases like "under $4/kg", "max $3.50", "budget $5".
- "mild", "skin", "sensitive", "conditioner", "shampoo", or "lotion" -> application_type=personal_care and min_perf >= 85.
- "industrial", "degreaser", "hard surface", "metal cleaner", or "solvent cleaner" -> application_type=industrial and min_perf >= 85.
- "agricultural", "crop", "adjuvant", "fertilizer", "herbicide", or "biostimulant" -> application_type=agricultural.
- "pharma", "pharmaceutical", "tablet", "capsule", "excipient", or "oral suspension" -> application_type=pharmaceutical.
- "food", "beverage", "GRAS", "bakery", or "flavor" -> application_type=food.
- "laundry", "detergent", "fabric", or "softener" -> application_type=fabric_laundry.
- "paint", "coating", "varnish", or "binder" -> application_type=paint_coatings.
- "98%", "99%", "100%", "fully bio", "bio-based", or "naturally derived" should raise min_bio appropriately.
- Be conservative and useful. If uncertain, choose the closest application_type and explain briefly in reasoning."""

_USER_PROMPT = 'Customer brief: """{text}"""'


def _log_error(message: str) -> None:
    print(message)
    try:
        with open("/tmp/intelliform_llm_errors.txt", "a") as handle:
            handle.write(message + "\n---\n")
    except Exception:
        pass


def _extract_json_object(text: str) -> str:
    cleaned = text.replace("```json", "").replace("```", "").strip()
    start = cleaned.find("{")
    end = cleaned.rfind("}")
    if start == -1 or end == -1 or end <= start:
        raise ValueError("No complete JSON object found in model output.")
    return cleaned[start:end + 1]


def _repair_json_text(text: str) -> str:
    repaired = _extract_json_object(text)
    repaired = repaired.replace("\u201c", '"').replace("\u201d", '"')
    repaired = repaired.replace("\u2018", "'").replace("\u2019", "'")
    repaired = re.sub(r",(\s*[}\]])", r"\1", repaired)
    return repaired


def _parse_json_payload(raw: str) -> dict:
    try:
        return json.loads(_extract_json_object(raw))
    except Exception:
        return json.loads(_repair_json_text(raw))


def _clamp(value: float, low: float, high: float) -> float:
    return max(low, min(high, value))


def _to_float(value, default: float) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _infer_application_type(text: str) -> str:
    t = text.lower()
    if any(term in t for term in ["conditioner", "shampoo", "lotion", "serum", "cleanser", "skin", "hair"]):
        return "personal_care"
    if any(term in t for term in ["laundry", "detergent", "fabric", "softener", "stain remover"]):
        return "fabric_laundry"
    if any(term in t for term in ["paint", "coating", "varnish", "primer"]):
        return "paint_coatings"
    if any(term in t for term in ["pharma", "pharmaceutical", "tablet", "capsule", "excipient", "oral suspension"]):
        return "pharmaceutical"
    if any(term in t for term in ["food", "beverage", "bakery", "flavor", "flavour", "gras"]):
        return "food"
    if any(term in t for term in ["agricultural", "crop", "adjuvant", "fertilizer", "biostimulant", "herbicide"]):
        return "agricultural"
    if any(term in t for term in ["industrial", "degreaser", "hard surface", "metal cleaner", "solvent cleaner"]):
        return "industrial"
    return "unknown"


def _normalize_application_type(value: str, raw_input: str) -> str:
    normalized = str(value or "").strip().lower().replace("-", "_")
    normalized = APPLICATION_TYPE_ALIASES.get(normalized, normalized)
    if normalized in VALID_APPLICATION_TYPES:
        return normalized
    inferred = _infer_application_type(raw_input)
    return inferred if inferred in VALID_APPLICATION_TYPES else "unknown"


def _build_result(data: dict, backend: str, raw_input: str) -> ParseResult:
    raw_text = raw_input.lower()
    application_type = _normalize_application_type(data.get("application_type", "unknown"), raw_input)
    profile = get_profile(application_type) if application_type != "unknown" else None

    default_cost = profile.default_max_cost if profile else 5.0
    default_bio = profile.default_min_bio if profile else 90.0
    default_perf = profile.default_min_perf if profile else 82.0

    high_bio_verticals = {"personal_care", "food"}
    medium_bio_verticals = {"fabric_laundry", "agricultural"}
    lower_bio_verticals = {"industrial", "paint_coatings"}

    if any(token in raw_text for token in ["100%", "99%", "98%", "fully bio"]):
        default_bio = max(default_bio, 98.0)
    elif any(token in raw_text for token in ["bio-based", "biobased", "naturally derived", "plant-based", "renewable"]):
        if application_type in high_bio_verticals:
            default_bio = max(default_bio, min(95.0, default_bio + 10.0))
        elif application_type in medium_bio_verticals:
            default_bio = max(default_bio, min(85.0, default_bio + 8.0))
        elif application_type in lower_bio_verticals:
            default_bio = max(default_bio, min(70.0, default_bio + 10.0))
        else:
            default_bio = max(default_bio, min(85.0, default_bio + 5.0))

    if any(token in raw_text for token in ["industrial", "degreaser", "cleaner", "mild", "skin", "conditioner", "shampoo"]):
        default_perf = max(default_perf, 85.0)

    if any(token in raw_text for token in ["release agent", "release coating", "concrete release", "demold"]):
        default_bio = min(default_bio, 70.0)
        default_perf = max(default_perf, 78.0)

    if any(token in raw_text for token in ["fabric softener", "softener"]):
        default_bio = min(max(default_bio, 55.0), 75.0)
        default_perf = max(default_perf, 72.0)

    if any(token in raw_text for token in ["detergent", "laundry liquid", "stain remover"]):
        default_bio = min(max(default_bio, 50.0), 80.0)
        default_perf = max(default_perf, 75.0)

    max_cost = _clamp(_to_float(data.get("max_cost"), default_cost), 0.5, 100.0)
    min_bio = _clamp(_to_float(data.get("min_bio"), default_bio), 0.0, 100.0)
    min_perf = _clamp(_to_float(data.get("min_perf"), default_perf), 0.0, 100.0)

    if application_type != "pharmaceutical":
        min_bio = max(min_bio, 80.0)
    min_perf = max(min_perf, 65.0)

    reasoning = str(data.get("reasoning", "")).strip()
    if not reasoning:
        reasoning = f"Parsed as {application_type.replace('_', ' ')} with cost, bio-based, and performance targets normalized for IntelliForm."

    return ParseResult(
        max_cost=round(max_cost, 2),
        min_bio=round(min_bio, 1),
        min_perf=round(min_perf, 1),
        application_type=application_type,
        reasoning=reasoning,
        parser_backend=backend,
        raw_input=raw_input,
    )


def _parse_with_groq(text: str) -> Optional[ParseResult]:
    api_key = os.getenv("GROQ_API_KEY", "")
    if not api_key:
        return None
    try:
        from groq import Groq

        client = Groq(api_key=api_key)
        model = os.getenv("GROQ_MODEL", "llama-3.3-70b-versatile")
        response = client.chat.completions.create(
            model=model,
            messages=[
                {"role": "system", "content": _SYSTEM_PROMPT},
                {"role": "user", "content": _USER_PROMPT.format(text=text)},
            ],
            temperature=0.1,
            max_tokens=300,
        )
        raw = response.choices[0].message.content or ""
        return _build_result(_parse_json_payload(raw), f"groq/{model}", text)
    except Exception as exc:
        _log_error(f"[llm_parser] Groq failed: {type(exc).__name__}: {exc}")
        return None


def _parse_with_anthropic(text: str) -> Optional[ParseResult]:
    api_key = os.getenv("ANTHROPIC_API_KEY", "")
    if not api_key:
        return None
    try:
        import anthropic

        client = anthropic.Anthropic(api_key=api_key)
        model = os.getenv("ANTHROPIC_MODEL", "claude-sonnet-4-5")
        response = client.messages.create(
            model=model,
            max_tokens=300,
            system=_SYSTEM_PROMPT,
            messages=[{"role": "user", "content": _USER_PROMPT.format(text=text)}],
        )
        raw = response.content[0].text.strip()
        return _build_result(_parse_json_payload(raw), f"anthropic/{model}", text)
    except Exception as exc:
        _log_error(f"[llm_parser] Anthropic failed: {type(exc).__name__}: {exc}\n{traceback.format_exc()}")
        return None


def _parse_with_openai(text: str) -> Optional[ParseResult]:
    api_key = os.getenv("OPENAI_API_KEY", "")
    if not api_key:
        return None
    try:
        from openai import OpenAI

        client = OpenAI(api_key=api_key)
        model = os.getenv("OPENAI_MODEL", "gpt-4o-mini")
        response = client.chat.completions.create(
            model=model,
            messages=[
                {"role": "system", "content": _SYSTEM_PROMPT},
                {"role": "user", "content": _USER_PROMPT.format(text=text)},
            ],
            temperature=0.1,
            max_tokens=300,
        )
        raw = response.choices[0].message.content or ""
        return _build_result(_parse_json_payload(raw), f"openai/{model}", text)
    except Exception as exc:
        _log_error(f"[llm_parser] OpenAI failed: {type(exc).__name__}: {exc}\n{traceback.format_exc()}")
        return None


def _parse_with_ollama(text: str) -> Optional[ParseResult]:
    host = os.getenv("OLLAMA_HOST", "http://localhost:11434")
    model = os.getenv("OLLAMA_MODEL", "llama3")
    try:
        resp = requests.post(
            f"{host}/api/generate",
            json={
                "model": model,
                "prompt": f"{_SYSTEM_PROMPT}\n\n{_USER_PROMPT.format(text=text)}",
                "stream": False,
                "format": "json",
            },
            timeout=20,
        )
        resp.raise_for_status()
        raw = resp.json().get("response", "")
        return _build_result(_parse_json_payload(raw), f"ollama/{model}", text)
    except Exception as exc:
        _log_error(f"[llm_parser] Ollama failed: {type(exc).__name__}: {exc}")
        return None


def _parse_with_regex(text: str) -> ParseResult:
    t = text.lower()
    m = re.search(r"(?:under|max|budget|below)\s*\$?\s*(\d+\.?\d*)\s*(?:/kg)?", t)
    app = _infer_application_type(text)
    profile = get_profile(app) if app != "unknown" else None
    max_cost = float(m.group(1)) if m else (profile.default_max_cost if profile else 5.0)
    base_bio = profile.default_min_bio if profile else 90.0
    base_perf = profile.default_min_perf if profile else 82.0

    if any(x in t for x in ["100%", "99%", "98%", "fully bio", "fully natural"]):
        min_bio = 98.0
    elif any(x in t for x in ["bio-based", "biobased", "naturally derived", "high bio", "sustainable", "plant-based"]):
        if app in {"personal_care", "food"}:
            min_bio = min(95.0, base_bio + 10.0)
        elif app in {"fabric_laundry", "agricultural"}:
            min_bio = min(85.0, base_bio + 8.0)
        elif app in {"industrial", "paint_coatings"}:
            min_bio = min(70.0, base_bio + 10.0)
        else:
            min_bio = min(85.0, base_bio + 5.0)
    else:
        min_bio = base_bio

    if any(x in t for x in ["release agent", "release coating", "concrete release", "demold"]):
        min_bio = min(min_bio, 70.0)
        base_perf = max(base_perf, 78.0)

    if any(x in t for x in ["fabric softener", "softener"]):
        min_bio = min(max(min_bio, 55.0), 75.0)
        base_perf = max(base_perf, 72.0)

    if any(x in t for x in ["detergent", "laundry liquid", "stain remover"]):
        min_bio = min(max(min_bio, 50.0), 80.0)
        base_perf = max(base_perf, 75.0)

    if any(x in t for x in ["foaming", "high foam", "surfactant"]):
        min_perf = 88.0
    elif any(x in t for x in ["skin", "cosmetics", "personal care", "mild", "industrial", "degreaser", "conditioner", "shampoo"]):
        min_perf = 85.0
    else:
        min_perf = base_perf

    return ParseResult(
        max_cost=round(max_cost, 2),
        min_bio=round(max(min_bio, 80.0 if app != "pharmaceutical" else min_bio), 1),
        min_perf=round(max(min_perf, 65.0), 1),
        application_type=app,
        reasoning="Extracted using regex fallback after local rule-based parsing.",
        parser_backend="regex",
        raw_input=text,
    )


def parse_request(text: str) -> ParseResult:
    """
    Parse a customer NL request into formulation constraints.

    Routing logic:
    - LLM_PROVIDER=groq|anthropic|openai|ollama -> use only that provider
    - LLM_PROVIDER=auto (default) -> try all in order, fall back to regex
    """
    provider = os.getenv("LLM_PROVIDER", "auto").lower()

    if provider == "groq":
        result = _parse_with_groq(text) or _parse_with_regex(text)
    elif provider == "anthropic":
        result = _parse_with_anthropic(text) or _parse_with_regex(text)
    elif provider == "openai":
        result = _parse_with_openai(text) or _parse_with_regex(text)
    elif provider == "ollama":
        result = _parse_with_ollama(text) or _parse_with_regex(text)
    else:
        result = (
            _parse_with_groq(text)
            or _parse_with_anthropic(text)
            or _parse_with_openai(text)
            or _parse_with_ollama(text)
            or _parse_with_regex(text)
        )

    track(
        "llm_parser_used",
        {
            "parser_backend": result.parser_backend,
            "application_type": result.application_type,
            "max_cost": result.max_cost,
            "min_bio": result.min_bio,
            "min_perf": result.min_perf,
            "input_length": len(text),
            "llm_provider_env": provider,
        },
    )

    return result
