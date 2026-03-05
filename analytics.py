"""
modules/analytics.py
PostHog analytics wrapper for IntelliForm v0.7.

All PostHog calls route through here so:
  - Analytics never crash the app (all try/except)
  - Event schema is centrally documented
  - Swapping providers later = one file change
"""
import os
import uuid
import streamlit as st

try:
    import posthog as _posthog
    _posthog.project_api_key = os.getenv(
        "POSTHOG_API_KEY",
        "phc_sCgIRYnbPtiYxwwxp4AodiswroSxNLPrVNvYOlpyq3z"   # Default project fallback
    )
    _posthog.host = os.getenv("POSTHOG_HOST", "https://us.i.posthog.com")
    POSTHOG_AVAILABLE = True
except ImportError:
    POSTHOG_AVAILABLE = False


# ── Session identity ─────────────────────────────────────────────────────────

def get_session_id() -> str:
    """
    Stable anonymous ID for this browser session.
    Lives in st.session_state so it survives Streamlit reruns
    but resets on a fresh browser tab.
    """
    if "session_id" not in st.session_state:
        st.session_state.session_id = str(uuid.uuid4())
    return st.session_state.session_id


def identify_user(email: str = None, name: str = None, company: str = None):
    """
    Optionally associate the session with a real user identity.
    Call this when the user voluntarily provides contact info
    (e.g. on the pilot booking form).

    PostHog will merge anonymous session events with this identity.
    """
    if not POSTHOG_AVAILABLE:
        return
    try:
        props = {}
        if name:    props["name"] = name
        if email:   props["email"] = email
        if company: props["company"] = company
        _posthog.identify(get_session_id(), props)
    except Exception:
        pass


# ── Core tracking function ───────────────────────────────────────────────────

def track(event: str, properties: dict = None):
    """
    Fire a PostHog event.

    EVENT CATALOGUE (keep this updated as you add events):
    ┌─────────────────────────────┬──────────────────────────────────────────┐
    │ Event                       │ When fired                               │
    ├─────────────────────────────┼──────────────────────────────────────────┤
    │ session_started             │ First load of the app                    │
    │ llm_parser_used             │ After LLM parse, records which backend   │
    │ formulation_generated       │ Successful PuLP optimization             │
    │ constraints_relaxed         │ First attempt infeasible, auto-relax     │
    │ error_optimization_failed   │ Infeasible even after relaxation         │
    │ pilot_button_clicked        │ User clicks Book NJ Pilot                │
    │ export_proposal             │ User downloads the proposal markdown     │
    └─────────────────────────────┴──────────────────────────────────────────┘
    """
    if not POSTHOG_AVAILABLE:
        return
    try:
        _posthog.capture(
            distinct_id=get_session_id(),
            event=event,
            properties=properties or {}
        )
    except Exception:
        pass  # Analytics must never crash the app
