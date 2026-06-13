"""
modules/api_key_store.py
Multi-tenant API key management for IntelliForm.

Keys are prefixed "if_" and 43 chars long (url-safe base64 of 32 random bytes).
Storage: /tmp/intelliform_api_keys.json (flat JSON; swap for Supabase table at scale).
Hashing:  SHA-256 of the raw key — raw key shown only once at creation.
Rate limiting: sliding 1-hour window, in-process counter (per-server, not distributed).

Supported scopes:
  feedback:write   — POST /api/v1/al/feedback, /api/v1/al/feedback/batch
  queue:read       — GET  /api/v1/al/experiment-queue, /api/v1/al/stats, /api/v1/al/feedback
  webhooks:write   — POST/DELETE /api/v1/webhooks
  sds:read         — POST /api/v1/export/sds
  qsar:read        — POST /api/v1/predict/qsar
  admin            — full access; only for ADMIN_API_KEY env override

Backward compat: if LIMS_API_KEY env var is set AND no key store entry matches,
the env var is tried as a master key with scope ["*"] (all scopes).
"""
from __future__ import annotations

import hashlib
import json
import os
import secrets
import time
from dataclasses import dataclass, asdict, field
from typing import Dict, List, Optional, Tuple

_STORE_PATH = "/tmp/intelliform_api_keys.json"
_RATE_COUNTERS: Dict[str, List[float]] = {}   # key_id → list of recent request timestamps

ALL_SCOPES = [
    "feedback:write", "queue:read", "webhooks:write",
    "sds:read", "qsar:read", "admin",
]


# ── Dataclass ─────────────────────────────────────────────────────────────────

@dataclass
class APIKey:
    id:                 str
    key_prefix:         str          # first 12 chars — safe to display
    key_hash:           str          # SHA-256 hex of the full raw key
    org_name:           str
    scopes:             List[str]    # subset of ALL_SCOPES, or ["*"] for master
    rate_limit_per_hour: int
    created_at:         str
    last_used_at:       Optional[str] = None
    active:             bool = True
    description:        str = ""


# ── Persistence ───────────────────────────────────────────────────────────────

def _load() -> List[APIKey]:
    if not os.path.exists(_STORE_PATH):
        return []
    try:
        with open(_STORE_PATH) as f:
            raw = json.load(f)
        return [APIKey(**r) for r in raw]
    except Exception:
        return []


def _save(keys: List[APIKey]) -> None:
    try:
        with open(_STORE_PATH, "w") as f:
            json.dump([asdict(k) for k in keys], f, indent=2)
    except Exception:
        pass


# ── Rate limiting ─────────────────────────────────────────────────────────────

def _within_rate_limit(key_id: str, limit_per_hour: int) -> bool:
    now          = time.time()
    window_start = now - 3600.0
    hits         = [t for t in _RATE_COUNTERS.get(key_id, []) if t > window_start]
    if len(hits) >= limit_per_hour:
        return False
    hits.append(now)
    _RATE_COUNTERS[key_id] = hits
    return True


# ── Public API ────────────────────────────────────────────────────────────────

def create_api_key(
    org_name:            str,
    scopes:              Optional[List[str]] = None,
    rate_limit_per_hour: int = 200,
    description:         str = "",
) -> Tuple[APIKey, str]:
    """
    Generate a new API key.
    Returns (APIKey, raw_key). raw_key is shown only once — store it securely.
    Raises ValueError for unknown scopes.
    """
    scopes = scopes or list(ALL_SCOPES)
    unknown = set(scopes) - set(ALL_SCOPES) - {"*"}
    if unknown:
        raise ValueError(f"Unknown scope(s): {unknown}. Valid: {ALL_SCOPES}")

    raw_key  = "if_" + secrets.token_urlsafe(32)
    key_hash = hashlib.sha256(raw_key.encode()).hexdigest()
    key_id   = secrets.token_hex(4)

    key = APIKey(
        id=key_id,
        key_prefix=raw_key[:12],
        key_hash=key_hash,
        org_name=org_name,
        scopes=scopes,
        rate_limit_per_hour=rate_limit_per_hour,
        created_at=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        active=True,
        description=description,
    )
    keys = _load()
    keys.append(key)
    _save(keys)
    return key, raw_key


def validate_api_key(raw_key: str) -> Optional[APIKey]:
    """
    Look up a raw key. Returns the APIKey if valid and active, else None.
    Updates last_used_at and enforces rate limit.
    Also checks LIMS_API_KEY env var as a master-key fallback.
    """
    if not raw_key:
        return None

    # Master key fallback (LIMS_API_KEY env var — backward compat)
    master = os.getenv("LIMS_API_KEY", "")
    if master and raw_key == master:
        return APIKey(
            id="__env_master__",
            key_prefix=raw_key[:8],
            key_hash="",
            org_name="env-master",
            scopes=["*"],
            rate_limit_per_hour=10_000,
            created_at="",
            active=True,
        )

    key_hash = hashlib.sha256(raw_key.encode()).hexdigest()
    keys     = _load()
    for key in keys:
        if key.key_hash == key_hash and key.active:
            if not _within_rate_limit(key.id, key.rate_limit_per_hour):
                return None   # rate limited — caller will raise 429
            # Update last_used_at
            key.last_used_at = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
            _save(keys)
            return key
    return None


def has_scope(key: APIKey, required_scope: str) -> bool:
    """Return True if key has the required scope (or wildcard '*')."""
    return "*" in key.scopes or "admin" in key.scopes or required_scope in key.scopes


def revoke_api_key(key_id: str) -> bool:
    """Deactivate a key by ID. Returns True if found."""
    keys   = _load()
    before = len([k for k in keys if k.active])
    for k in keys:
        if k.id == key_id:
            k.active = False
    _save(keys)
    after = len([k for k in keys if k.active])
    return after < before


def list_api_keys() -> List[APIKey]:
    """Return all registered keys (active and inactive)."""
    return _load()


def get_api_key(key_id: str) -> Optional[APIKey]:
    return next((k for k in _load() if k.id == key_id), None)
