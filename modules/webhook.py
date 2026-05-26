"""
modules/webhook.py
Webhook egress for IntelliForm — fires POST events to registered LIMS/CI URLs
after model retrains, experiment queue updates, or batch feedback ingestion.

Storage: /tmp/intelliform_webhooks.json (flat JSON; Supabase upgrade path available).
Security: each delivery is signed with HMAC-SHA256(payload, LIMS_API_KEY).
Retry:    3 attempts with 1s / 2s / 4s back-off. Failures are logged, not raised.

Supported events:
  model.retrained        — QSAR GBR retrained after AL feedback
  stability.retrained    — Stability GBR retrained after lab feedback
  batch.ingested         — Batch feedback endpoint accepted records
  test.ping              — Fired by POST /api/v1/webhooks/test

Payload schema (all events):
  {
    "event":         str,
    "timestamp":     ISO-8601 str,
    "al_round":      int | None,
    "total_feedback": int | None,
    "source":        str,           # "qsar" | "stability" | "api"
    "details":       dict           # event-specific extra fields
  }
"""
from __future__ import annotations

import hashlib
import hmac
import json
import os
import time
import uuid
from dataclasses import dataclass, asdict, field
from typing import Dict, List, Optional

import requests

_STORE_PATH   = "/tmp/intelliform_webhooks.json"
_REQUEST_TIMEOUT = 6      # seconds per attempt
_MAX_RETRIES     = 3
_BACKOFF_BASE    = 1.0    # seconds; doubled each retry

_SUPPORTED_EVENTS = {
    "model.retrained",
    "stability.retrained",
    "batch.ingested",
    "test.ping",
}


# ── Data schema ───────────────────────────────────────────────────────────────

@dataclass
class WebhookConfig:
    id:         str
    url:        str
    events:     List[str]          # list of event names this hook subscribes to
    created_at: str
    active:     bool = True
    description: str = ""


# ── Persistence ───────────────────────────────────────────────────────────────

def _load() -> List[WebhookConfig]:
    if not os.path.exists(_STORE_PATH):
        return []
    try:
        with open(_STORE_PATH) as f:
            raw = json.load(f)
        return [WebhookConfig(**r) for r in raw]
    except Exception:
        return []


def _save(hooks: List[WebhookConfig]) -> None:
    try:
        with open(_STORE_PATH, "w") as f:
            json.dump([asdict(h) for h in hooks], f, indent=2)
    except Exception:
        pass


# ── HMAC signing ──────────────────────────────────────────────────────────────

def _sign(payload_bytes: bytes) -> str:
    secret = os.getenv("LIMS_API_KEY", "intelliform-dev-secret").encode()
    return "sha256=" + hmac.new(secret, payload_bytes, hashlib.sha256).hexdigest()


# ── Delivery ──────────────────────────────────────────────────────────────────

def _deliver(hook: WebhookConfig, payload: dict) -> bool:
    body = json.dumps(payload, default=str).encode()
    sig  = _sign(body)
    headers = {
        "Content-Type":           "application/json",
        "X-IntelliForm-Event":    payload.get("event", ""),
        "X-IntelliForm-Delivery": str(uuid.uuid4()),
        "X-IntelliForm-Signature": sig,
    }
    delay = _BACKOFF_BASE
    for attempt in range(_MAX_RETRIES):
        try:
            resp = requests.post(
                hook.url, data=body, headers=headers,
                timeout=_REQUEST_TIMEOUT, allow_redirects=False,
            )
            if resp.status_code < 400:
                return True
            print(f"[webhook] {hook.url} → HTTP {resp.status_code} (attempt {attempt+1})")
        except Exception as exc:
            print(f"[webhook] {hook.url} delivery error (attempt {attempt+1}): {exc}")
        if attempt < _MAX_RETRIES - 1:
            time.sleep(delay)
            delay *= 2
    return False


# ── Public API ────────────────────────────────────────────────────────────────

def register_webhook(
    url:         str,
    events:      Optional[List[str]] = None,
    description: str = "",
) -> WebhookConfig:
    """
    Register a new webhook URL.
    events defaults to all supported events if not specified.
    Raises ValueError for invalid URLs or event names.
    """
    if not url.startswith(("http://", "https://")):
        raise ValueError(f"Webhook URL must start with http:// or https://: {url!r}")
    events = events or list(_SUPPORTED_EVENTS)
    unknown = set(events) - _SUPPORTED_EVENTS
    if unknown:
        raise ValueError(f"Unknown event(s): {unknown}. Supported: {_SUPPORTED_EVENTS}")

    hook = WebhookConfig(
        id=str(uuid.uuid4())[:8],
        url=url,
        events=events,
        created_at=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        active=True,
        description=description,
    )
    hooks = _load()
    hooks.append(hook)
    _save(hooks)
    return hook


def delete_webhook(webhook_id: str) -> bool:
    """Remove a webhook by ID. Returns True if found and removed."""
    hooks = _load()
    before = len(hooks)
    hooks  = [h for h in hooks if h.id != webhook_id]
    _save(hooks)
    return len(hooks) < before


def list_webhooks() -> List[WebhookConfig]:
    """Return all registered webhook configs."""
    return _load()


def fire_event(
    event:   str,
    source:  str = "api",
    details: Optional[dict] = None,
    al_round:       Optional[int] = None,
    total_feedback: Optional[int] = None,
) -> Dict[str, int]:
    """
    Fire an event to all subscribed, active webhooks.
    Returns {"fired": N, "failed": M}.
    Non-blocking best-effort: failures are printed, never raised.
    """
    if event not in _SUPPORTED_EVENTS:
        print(f"[webhook] Unknown event '{event}' — skipped")
        return {"fired": 0, "failed": 0}

    payload = {
        "event":          event,
        "timestamp":      time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "source":         source,
        "al_round":       al_round,
        "total_feedback": total_feedback,
        "details":        details or {},
    }

    hooks  = [h for h in _load() if h.active and event in h.events]
    fired  = failed = 0
    for hook in hooks:
        ok = _deliver(hook, payload)
        if ok:
            fired += 1
        else:
            failed += 1

    if hooks:
        print(f"[webhook] event={event} fired={fired} failed={failed}")
    return {"fired": fired, "failed": failed}
