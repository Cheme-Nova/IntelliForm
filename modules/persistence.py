"""
modules/persistence.py
Supabase persistence layer for IntelliForm v0.9.

Provides project storage, retrieval, and active learning feedback
across sessions. Falls back gracefully to in-memory storage when
Supabase credentials are not configured.

Supabase schema (run migrations/001_create_tables.sql to set up):

  Table: intelliform_projects
    id              uuid  PRIMARY KEY DEFAULT uuid_generate_v4()
    session_id      text  NOT NULL
    user_email      text
    created_at      timestamptz DEFAULT now()
    application     text
    blend           jsonb
    cost_per_kg     float
    bio_pct         float
    perf_score      float
    eco_score       float
    eco_grade       text
    optimizer       text
    parser          text
    relaxed         boolean
    nl_input        text

  Table: intelliform_feedback
    id              uuid  PRIMARY KEY DEFAULT uuid_generate_v4()
    session_id      text
    smiles          text
    target          text   -- 'Biodegradability' | 'Ecotoxicity' | 'Performance'
    predicted       float
    actual          float
    created_at      timestamptz DEFAULT now()

  Table: intelliform_pilot_bookings
    id              uuid  PRIMARY KEY DEFAULT uuid_generate_v4()
    session_id      text
    user_email      text
    blend           jsonb
    cost_per_kg     float
    batch_kg        int
    quote_usd       float
    application     text
    created_at      timestamptz DEFAULT now()
    status          text DEFAULT 'pending'
"""
import os
import json
from datetime import datetime
from typing import Dict, List, Optional, Any

# ── Optional Supabase import ──────────────────────────────────────────────────
try:
    from supabase import create_client, Client
    SUPABASE_OK = True
except ImportError:
    SUPABASE_OK = False

# ── In-memory fallback store ──────────────────────────────────────────────────
_MEMORY_STORE: Dict[str, List] = {
    "projects": [],
    "feedback": [],
    "bookings": [],
}


# ── Client singleton ──────────────────────────────────────────────────────────

def _get_client() -> Optional[Any]:
    if not SUPABASE_OK:
        return None
    url  = os.getenv("SUPABASE_URL", "")
    key  = os.getenv("SUPABASE_ANON_KEY", "")
    if not url or not key:
        return None
    try:
        return create_client(url, key)
    except Exception:
        return None


def is_connected() -> bool:
    """Returns True if Supabase is configured and reachable."""
    return _get_client() is not None


# ── Project storage ───────────────────────────────────────────────────────────

def save_project(project: Dict, session_id: str) -> bool:
    """
    Persist a formulation project.
    Returns True on success, False on error (never raises).
    """
    record = {
        "session_id":  session_id,
        "created_at":  datetime.utcnow().isoformat(),
        "application": project.get("application"),
        "blend":       project.get("blend", {}),
        "cost_per_kg": project.get("cost"),
        "bio_pct":     project.get("bio"),
        "perf_score":  project.get("perf"),
        "eco_score":   project.get("eco_score"),
        "eco_grade":   project.get("eco_grade"),
        "optimizer":   project.get("optimizer", "pulp"),
        "parser":      project.get("parser", "regex"),
        "relaxed":     project.get("relaxed", False),
        "nl_input":    project.get("input", ""),
    }

    client = _get_client()
    if client:
        try:
            client.table("intelliform_projects").insert(record).execute()
            return True
        except Exception as e:
            print(f"[persistence] Supabase save failed: {e}")

    # Fallback to memory
    _MEMORY_STORE["projects"].append(record)
    return False  # indicates memory fallback


def load_projects(session_id: str, limit: int = 50) -> List[Dict]:
    """
    Load projects for a session (most recent first).
    """
    client = _get_client()
    if client:
        try:
            resp = (
                client.table("intelliform_projects")
                .select("*")
                .eq("session_id", session_id)
                .order("created_at", desc=True)
                .limit(limit)
                .execute()
            )
            return resp.data or []
        except Exception as e:
            print(f"[persistence] Supabase load failed: {e}")

    # Memory fallback
    return [p for p in reversed(_MEMORY_STORE["projects"])
            if p.get("session_id") == session_id][:limit]


def load_all_projects(limit: int = 200) -> List[Dict]:
    """Load all projects across sessions (admin view)."""
    client = _get_client()
    if client:
        try:
            resp = (
                client.table("intelliform_projects")
                .select("*")
                .order("created_at", desc=True)
                .limit(limit)
                .execute()
            )
            return resp.data or []
        except Exception as e:
            print(f"[persistence] Supabase load_all failed: {e}")

    return list(reversed(_MEMORY_STORE["projects"]))[:limit]


# ── Active learning feedback ──────────────────────────────────────────────────

def save_feedback(session_id: str, smiles: str, target: str,
                  predicted: float, actual: float) -> bool:
    """Store a validated data point for QSAR active learning."""
    record = {
        "session_id": session_id,
        "smiles":     smiles,
        "target":     target,
        "predicted":  predicted,
        "actual":     actual,
        "created_at": datetime.utcnow().isoformat(),
    }
    client = _get_client()
    if client:
        try:
            client.table("intelliform_feedback").insert(record).execute()
            return True
        except Exception as e:
            print(f"[persistence] Feedback save failed: {e}")

    _MEMORY_STORE["feedback"].append(record)
    return False


def load_feedback_for_retraining() -> List[Dict]:
    """Load all validated feedback for model retraining."""
    client = _get_client()
    if client:
        try:
            resp = client.table("intelliform_feedback").select("*").execute()
            return resp.data or []
        except Exception:
            pass
    return _MEMORY_STORE["feedback"]


# ── Pilot bookings ────────────────────────────────────────────────────────────

def save_booking(session_id: str, blend: Dict, cost_per_kg: float,
                 batch_kg: int, quote_usd: float, application: str,
                 user_email: str = "") -> bool:
    """Record a pilot booking request."""
    record = {
        "session_id":  session_id,
        "user_email":  user_email,
        "blend":       blend,
        "cost_per_kg": cost_per_kg,
        "batch_kg":    batch_kg,
        "quote_usd":   quote_usd,
        "application": application,
        "created_at":  datetime.utcnow().isoformat(),
        "status":      "pending",
    }
    client = _get_client()
    if client:
        try:
            client.table("intelliform_pilot_bookings").insert(record).execute()
            return True
        except Exception as e:
            print(f"[persistence] Booking save failed: {e}")

    _MEMORY_STORE["bookings"].append(record)
    return False


# ── SQL migration helper ──────────────────────────────────────────────────────

MIGRATION_SQL = """
-- Run this in your Supabase SQL editor to set up IntelliForm tables
-- migrations/001_create_tables.sql

CREATE EXTENSION IF NOT EXISTS "uuid-ossp";

CREATE TABLE IF NOT EXISTS intelliform_projects (
    id          UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    session_id  TEXT NOT NULL,
    user_email  TEXT,
    created_at  TIMESTAMPTZ DEFAULT NOW(),
    application TEXT,
    blend       JSONB,
    cost_per_kg FLOAT,
    bio_pct     FLOAT,
    perf_score  FLOAT,
    eco_score   FLOAT,
    eco_grade   TEXT,
    optimizer   TEXT,
    parser      TEXT,
    relaxed     BOOLEAN,
    nl_input    TEXT
);

CREATE TABLE IF NOT EXISTS intelliform_feedback (
    id          UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    session_id  TEXT,
    smiles      TEXT,
    target      TEXT,
    predicted   FLOAT,
    actual      FLOAT,
    created_at  TIMESTAMPTZ DEFAULT NOW()
);

CREATE TABLE IF NOT EXISTS intelliform_pilot_bookings (
    id          UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    session_id  TEXT,
    user_email  TEXT,
    blend       JSONB,
    cost_per_kg FLOAT,
    batch_kg    INT,
    quote_usd   FLOAT,
    application TEXT,
    created_at  TIMESTAMPTZ DEFAULT NOW(),
    status      TEXT DEFAULT 'pending'
);

-- Index for fast session lookups
CREATE INDEX IF NOT EXISTS idx_projects_session  ON intelliform_projects(session_id);
CREATE INDEX IF NOT EXISTS idx_feedback_session  ON intelliform_feedback(session_id);
CREATE INDEX IF NOT EXISTS idx_bookings_session  ON intelliform_pilot_bookings(session_id);
"""
