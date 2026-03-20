"""
modules/tiers.py
IntelliForm Tier System — honest version

Gate ONLY on things that actually cost ChemeNova money:
  - Anthropic / OpenAI API calls (per-token cost)
  - Custom DB upload (storage cost, support burden)
  - API access (infrastructure cost)
  - SSO (implementation cost)
  - Dedicated support (time cost)

Everything that costs $0 to provide is FREE on all tiers:
  - All optimization modes (LP, Pareto, Bayesian) — free libraries
  - All 7 verticals + constraints — code only
  - Pharma deep dive — code only
  - Mordred 1613 descriptors — free library
  - PubChemPy enrichment — free API
  - Carbon credits, EcoMetrics, stability — code only
  - Regulatory reports — code only
  - PDF export — reportlab is free
  - Groq LLM — free tier
  - PostHog analytics — free tier
  - Supabase persistence — free tier
"""
from __future__ import annotations
import os
from dataclasses import dataclass
from typing import List, Optional


@dataclass
class TierConfig:
    name: str
    display_name: str
    color: str
    emoji: str
    price: str
    # What actually costs money to provide
    llms_allowed: List[str]       # Anthropic/OpenAI cost per token
    max_batch_kg: int             # affects proposal PDF complexity only
    white_label: bool             # custom branding = support time
    api_access: bool              # needs infrastructure
    custom_db_upload: bool        # storage + validation support
    sso_enabled: bool             # implementation cost
    dedicated_support: bool       # time cost
    export_formats: List[str]
    features_summary: List[str]


TIERS = {
    "free": TierConfig(
        name="free",
        display_name="Free",
        color="#0D9488",
        emoji="🌱",
        price="$0 · MIT License · Deploy anywhere",
        llms_allowed=["groq"],          # Groq is free
        max_batch_kg=5000,
        white_label=False,
        api_access=False,
        custom_db_upload=False,
        sso_enabled=False,
        dedicated_support=False,
        export_formats=["csv", "pdf"],
        features_summary=[
            "✅ 1,197+ ingredient database",
            "✅ LP + Pareto + Bayesian optimization",
            "✅ All 7 verticals + specific constraints",
            "✅ Pharma deep dive (ICH/BCS/USP/NF)",
            "✅ QSAR ML — Mordred 1613 descriptors",
            "✅ PubChem auto-enrichment",
            "✅ Regulatory intelligence (7 frameworks)",
            "✅ Carbon credit calculator",
            "✅ EcoMetrics + stability prediction",
            "✅ PDF + CSV export",
            "✅ Groq LLM — free inference",
            "✅ Supabase persistence",
            "✅ MIT license — fork and deploy",
        ],
    ),

    "pro": TierConfig(
        name="pro",
        display_name="Pro",
        color="#D97706",
        emoji="⚡",
        price="$99 / month",
        llms_allowed=["groq", "anthropic", "openai"],  # these cost money
        max_batch_kg=50000,
        white_label=True,
        api_access=False,
        custom_db_upload=True,          # needs validation + storage
        sso_enabled=False,
        dedicated_support=False,
        export_formats=["csv", "pdf", "json", "excel"],
        features_summary=[
            "✅ Everything in Free",
            "✅ Claude 3.5 Sonnet + GPT-4o LLMs",
            "✅ Custom ingredient DB upload",
            "✅ White-label PDF branding",
            "✅ Excel + JSON export",
            "✅ 50,000 kg batch size",
            "✅ Priority email support",
        ],
    ),

    "enterprise": TierConfig(
        name="enterprise",
        display_name="Enterprise",
        color="#0A1628",
        emoji="🏢",
        price="Custom — contact shehan@chemenova.com",
        llms_allowed=["groq", "anthropic", "openai", "azure", "custom"],
        max_batch_kg=10_000_000,
        white_label=True,
        api_access=True,                # needs infra
        custom_db_upload=True,
        sso_enabled=True,               # needs implementation
        dedicated_support=True,         # needs time
        export_formats=["csv", "pdf", "json", "excel", "xml", "edi"],
        features_summary=[
            "✅ Everything in Pro",
            "✅ REST API access",
            "✅ SSO / SAML authentication",
            "✅ Dedicated cloud infrastructure",
            "✅ ERP / LIMS / SAP connectors",
            "✅ On-premise deployment",
            "✅ SLA 99.9% uptime",
            "✅ Dedicated Customer Success Manager",
            "✅ Custom vertical development",
            "✅ IP indemnification",
            "✅ Audit log (21 CFR Part 11)",
            "✅ Unlimited seats",
        ],
    ),
}


def get_tier(name: Optional[str] = None) -> TierConfig:
    if name is None:
        name = os.environ.get("INTELLIFORM_TIER", "free").lower()
    return TIERS.get(name, TIERS["free"])


def gate(tier: TierConfig, feature: str, label: str = "") -> bool:
    """
    Returns True if feature available on this tier.
    Shows upgrade prompt if not.
    Only gates features that actually cost money to provide.
    """
    ok = getattr(tier, feature, True)  # default True — don't gate by default
    if not ok and label:
        import streamlit as st
        next_t = _next_tier(tier)
        st.warning(
            f"🔒 **{label}** requires **{next_t}**. "
            f"[Upgrade →](mailto:shehan@chemenova.com"
            f"?subject=IntelliForm%20{next_t}%20Upgrade)"
        )
    return ok


def _next_tier(tier: TierConfig) -> str:
    order = ["free", "pro", "enterprise"]
    idx = order.index(tier.name) if tier.name in order else 0
    return TIERS[order[min(idx + 1, len(order) - 1)]].display_name
