"""
modules/analytics.py
PostHog analytics wrapper for IntelliForm.
"""
import os
import uuid
import streamlit as st

try:
    import posthog as _posthog
    _posthog.project_api_key = os.getenv("POSTHOG_API_KEY", "phc_sCgIRYnbPtiYxwwxp4AodiswroSxNLPrVNvYOlpyq3z")
    _posthog.host = os.getenv("POSTHOG_HOST", "https://us.i.posthog.com")
    POSTHOG_AVAILABLE = True
except ImportError:
    POSTHOG_AVAILABLE = False


def get_session_id() -> str:
    if "session_id" not in st.session_state:
        st.session_state.session_id = str(uuid.uuid4())
    return st.session_state.session_id


def identify_user(email: str = None, name: str = None, company: str = None):
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


def track(event: str, properties: dict = None):
    if not POSTHOG_AVAILABLE:
        return
    try:
        _posthog.capture(distinct_id=get_session_id(), event=event, properties=properties or {})
    except Exception:
        pass
