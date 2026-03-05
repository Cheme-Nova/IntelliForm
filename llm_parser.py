"""
modules/llm_parser.py
Natural language → formulation constraints parser for IntelliForm v0.7.

Priority chain:
  1. Groq  (llama-3.1-8b-instant, ~300ms, free tier generous)
  2. Ollama local (llama3 or mistral, zero latency after first load)
  3. Regex fallback (v0.6 behaviour, always works)

The parser returns a ParseResult dataclass so callers get typed fields
and the reasoning the LLM used — we surface this in the UI.
"""
import os
import re
import json
import requests
from dataclasses import dataclass, field
from typing import Optional

from modules.analytics import track


# ── Result schema ─────────────────────────────────────────────────────────────

@dataclass
class ParseResult:
    max_cost: float          # USD/kg ceiling
    min_bio: float           # bio-based % floor
    min_perf: float          # performance score floor
    application_type: str    # e.g. "cosmetics", "industrial", "food-safe"
    reasoning: str           # LLM chain-of-thought or "regex fallback"
    parser_backend: str      # "groq" | "ollama" | "regex"
    raw_input: str           # original user text


# ── Shared LLM system prompt ──────────────────────────────────────────────────

_SYSTEM_PROMPT = """You are a green chemistry formulation assistant for IntelliForm, an AI platform \
used by chemical engineers to design sustainable surfactant blends.

Your job: extract formulation constraints from a customer request.

Respond ONLY with a valid JSON object — no markdown, no preamble, no trailing text.

JSON schema:
{
  "max_cost":        <float, USD per kg, default 5.0>,
  "min_bio":         <float, bio-based percentage 0-100, default 95>,
  "min_perf":        <float, performance score 0-100, default 82>,
  "application_type":<string, one of: cosmetics | personal_care | industrial | food_safe | textile | unknown>,
  "reasoning":       <string, ≤2 sentences explaining your extraction>
}

Rules:
- If the user mentions "mild" or "skin", set min_perf ≥ 85
- If the user mentions "foaming" or "surfactant", set min_perf ≥ 88  
- If the user mentions "98%", "99%", "100%", "fully bio" or "fully natural", set min_bio ≥ 98
- Extract cost ceilings from phrases like "under $4/kg", "max $3.50", "budget $5"
- If no cost mentioned, default max_cost to 5.0
- Never set min_bio below 80 or min_perf below 65
- application_type should reflect the end-use even if not stated explicitly
"""

_USER_PROMPT_TEMPLATE = 'Customer request: """{text}"""'


# ── Groq backend ──────────────────────────────────────────────────────────────

def _parse_with_groq(text: str) -> Optional[ParseResult]:
    api_key = os.getenv("GROQ_API_KEY", "")
    if not api_key:
        return None
    try:
        from groq import Groq
        client = Groq(api_key=api_key)
        response = client.chat.completions.create(
            model="llama-3.1-8b-instant",
            messages=[
                {"role": "system", "content": _SYSTEM_PROMPT},
                {"role": "user",   "content": _USER_PROMPT_TEMPLATE.format(text=text)}
            ],
            temperature=0.1,      # low temp = deterministic extraction
            max_tokens=300,
            response_format={"type": "json_object"}
        )
        raw = response.choices[0].message.content
        data = json.loads(raw)
        return ParseResult(
            max_cost        = float(data.get("max_cost", 5.0)),
            min_bio         = float(data.get("min_bio", 95)),
            min_perf        = float(data.get("min_perf", 82)),
            application_type= str(data.get("application_type", "unknown")),
            reasoning       = str(data.get("reasoning", "")),
            parser_backend  = "groq",
            raw_input       = text
        )
    except Exception as e:
        # Log but don't crash — fall through to next backend
        print(f"[llm_parser] Groq failed: {e}")
        return None


# ── Ollama backend ────────────────────────────────────────────────────────────

def _parse_with_ollama(text: str) -> Optional[ParseResult]:
    host = os.getenv("OLLAMA_HOST", "http://localhost:11434")
    model = os.getenv("OLLAMA_MODEL", "llama3")
    try:
        payload = {
            "model": model,
            "prompt": f"{_SYSTEM_PROMPT}\n\n{_USER_PROMPT_TEMPLATE.format(text=text)}",
            "stream": False,
            "format": "json"
        }
        resp = requests.post(f"{host}/api/generate", json=payload, timeout=15)
        resp.raise_for_status()
        raw = resp.json().get("response", "")
        data = json.loads(raw)
        return ParseResult(
            max_cost        = float(data.get("max_cost", 5.0)),
            min_bio         = float(data.get("min_bio", 95)),
            min_perf        = float(data.get("min_perf", 82)),
            application_type= str(data.get("application_type", "unknown")),
            reasoning       = str(data.get("reasoning", "")),
            parser_backend  = "ollama",
            raw_input       = text
        )
    except Exception as e:
        print(f"[llm_parser] Ollama failed: {e}")
        return None


# ── Regex fallback (v0.6 behaviour, always succeeds) ─────────────────────────

def _parse_with_regex(text: str) -> ParseResult:
    t = text.lower()

    # Cost ceiling
    m = re.search(r'(?:under|max|budget|below)\s*\$?\s*(\d+\.?\d*)\s*(?:/kg)?', t)
    max_cost = float(m.group(1)) if m else 5.0

    # Bio-based floor
    if any(x in t for x in ['100%', '99%', '98%', 'fully bio', 'fully natural', 'high bio']):
        min_bio = 98.0
    else:
        min_bio = 95.0

    # Performance floor
    if any(x in t for x in ['foaming', 'high foam', 'surfactant']):
        min_perf = 88.0
    elif any(x in t for x in ['skin', 'cosmetics', 'personal care', 'mild']):
        min_perf = 85.0
    else:
        min_perf = 82.0

    # Application type
    if any(x in t for x in ['cosmetics', 'skin', 'personal care', 'lotion', 'cream']):
        app = "cosmetics"
    elif any(x in t for x in ['food', 'edible', 'fda']):
        app = "food_safe"
    elif any(x in t for x in ['textile', 'fabric', 'laundry']):
        app = "textile"
    elif any(x in t for x in ['industrial', 'cleaning', 'degreaser']):
        app = "industrial"
    else:
        app = "unknown"

    return ParseResult(
        max_cost        = max_cost,
        min_bio         = min_bio,
        min_perf        = min_perf,
        application_type= app,
        reasoning       = "Extracted using regex fallback (no LLM API key configured).",
        parser_backend  = "regex",
        raw_input       = text
    )


# ── Public interface ──────────────────────────────────────────────────────────

def parse_request(text: str) -> ParseResult:
    """
    Parse a customer natural-language request into formulation constraints.
    Tries Groq → Ollama → regex, in that order.
    Always succeeds (regex is guaranteed fallback).
    Fires PostHog event 'llm_parser_used' with backend info.
    """
    result = _parse_with_groq(text) or _parse_with_ollama(text) or _parse_with_regex(text)

    track("llm_parser_used", {
        "parser_backend":   result.parser_backend,
        "application_type": result.application_type,
        "max_cost":         result.max_cost,
        "min_bio":          result.min_bio,
        "min_perf":         result.min_perf,
        "input_length":     len(text)
    })

    return result
