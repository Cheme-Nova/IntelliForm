"""
modules/tiers.py
Tier enforcement for IntelliForm v1.0.

Tiers:
  free        — 3 runs/session, basic features, watermarked output
  pro         — $299/month, unlimited, all features
  team        — $899/month, 5 seats, custom ingredients, priority booking
  enterprise  — custom, white-label, API access

Set INTELLIFORM_TIER=free|pro|team|enterprise in Streamlit secrets.
Default is 'free' if not set.
"""
import os
import streamlit as st

# ── Constants ─────────────────────────────────────────────────────────────────
FREE_RUN_LIMIT = 3

TIER_CONFIG = {
    "free": {
        "name": "Free",
        "color": "#64748b",
        "price": "$0",
        "run_limit": FREE_RUN_LIMIT,
        "features": {
            "pareto":         False,
            "pdf":            False,
            "carbon":         False,
            "stability":      False,
            "comparison":     False,
            "llm_switch":     False,
            "pilot_booking":  False,
            "multi_llm":      False,
        }
    },
    "pro": {
        "name": "Professional",
        "color": "#0D9488",
        "price": "$299/mo",
        "run_limit": None,  # unlimited
        "features": {
            "pareto":         True,
            "pdf":            True,
            "carbon":         True,
            "stability":      True,
            "comparison":     True,
            "llm_switch":     True,
            "pilot_booking":  True,
            "multi_llm":      True,
        }
    },
    "team": {
        "name": "Team",
        "color": "#D97706",
        "price": "$899/mo",
        "run_limit": None,
        "features": {
            "pareto":         True,
            "pdf":            True,
            "carbon":         True,
            "stability":      True,
            "comparison":     True,
            "llm_switch":     True,
            "pilot_booking":  True,
            "multi_llm":      True,
            "custom_ingredients": True,
            "priority_booking":   True,
        }
    },
    "enterprise": {
        "name": "Enterprise",
        "color": "#0A1628",
        "price": "Custom",
        "run_limit": None,
        "features": {
            "pareto":         True,
            "pdf":            True,
            "carbon":         True,
            "stability":      True,
            "comparison":     True,
            "llm_switch":     True,
            "pilot_booking":  True,
            "multi_llm":      True,
            "custom_ingredients": True,
            "priority_booking":   True,
            "white_label":        True,
            "api_access":         True,
        }
    }
}

UPGRADE_URL = "mailto:shehan@chemenova.com?subject=IntelliForm Pro Upgrade"


# ── Core tier helpers ─────────────────────────────────────────────────────────

def get_tier() -> str:
    """Get current tier from environment. Default free."""
    return os.getenv("INTELLIFORM_TIER", "free").lower()


def is_free() -> bool:
    return get_tier() == "free"


def is_pro() -> bool:
    return get_tier() in ["pro", "team", "enterprise"]


def has_feature(feature: str) -> bool:
    tier = get_tier()
    config = TIER_CONFIG.get(tier, TIER_CONFIG["free"])
    return config["features"].get(feature, False)


# ── Run count management ──────────────────────────────────────────────────────

def get_run_count() -> int:
    return st.session_state.get("run_count", 0)


def increment_run_count():
    st.session_state["run_count"] = get_run_count() + 1


def runs_remaining() -> int:
    if not is_free():
        return 999
    return max(0, FREE_RUN_LIMIT - get_run_count())


def can_run() -> bool:
    if not is_free():
        return True
    return get_run_count() < FREE_RUN_LIMIT


# ── Gate functions (return True if allowed) ───────────────────────────────────

def gate_pareto() -> bool:       return has_feature("pareto")
def gate_pdf() -> bool:          return has_feature("pdf")
def gate_carbon() -> bool:       return has_feature("carbon")
def gate_stability() -> bool:    return has_feature("stability")
def gate_comparison() -> bool:   return has_feature("comparison")
def gate_llm_switch() -> bool:   return has_feature("llm_switch")
def gate_pilot_booking() -> bool:return has_feature("pilot_booking")


# ── UI Components ─────────────────────────────────────────────────────────────

def tier_badge() -> str:
    tier = get_tier()
    config = TIER_CONFIG.get(tier, TIER_CONFIG["free"])
    return f"{config['name']} ({config['price']})"


def show_upgrade_banner():
    """Show upgrade prompt when free user hits run limit."""
    st.error(f"""
**You've used all {FREE_RUN_LIMIT} free formulation runs for this session.**

Upgrade to **IntelliForm Professional** for unlimited runs, PDF proposals,
Pareto optimization, stability prediction, carbon credits, and pilot batch booking.

📧 [Contact us to upgrade]({UPGRADE_URL}) · $299/month · Cancel anytime
""")


def show_watermark():
    """Show free tier watermark."""
    if is_free():
        st.caption(
            "🆓 Free tier — [Upgrade to Pro]({}) for PDF export, Pareto, "
            "stability prediction, carbon credits & more".format(UPGRADE_URL)
        )


def show_locked_feature(feature_name: str, description: str = ""):
    """Show locked feature placeholder with upgrade prompt."""
    st.info(
        f"🔒 **{feature_name}** is available on IntelliForm Professional ($299/mo)\n\n"
        f"{description}\n\n"
        f"[Contact us to upgrade]({UPGRADE_URL})"
    )


def show_run_counter():
    """Show remaining runs for free tier."""
    if is_free():
        remaining = runs_remaining()
        if remaining > 0:
            st.caption(f"🆓 Free tier · {remaining} run{'s' if remaining != 1 else ''} remaining this session")
        else:
            show_upgrade_banner()
            st.stop()
