"""
modules/supply_chain.py
Supplier marketplace + supply-chain risk intelligence for IntelliForm.

Supplier keys:  ifs_ prefix (44 chars, url-safe base64 of 32 bytes)
Storage:        /tmp/intelliform_suppliers.json
                /tmp/intelliform_listings.json
Auth header:    X-Supplier-Key on supplier-facing endpoints

Supplier lifecycle:  register (pending) → admin approve (active) → gets ifs_ key
"""
from __future__ import annotations

import hashlib
import json
import os
import secrets
import time
from dataclasses import dataclass, asdict
from typing import Dict, List, Optional, Tuple

_SUPPLIERS_PATH = "/tmp/intelliform_suppliers.json"
_LISTINGS_PATH  = "/tmp/intelliform_listings.json"

SUPPORTED_CERTS = [
    "ISO9001", "ISO14001", "GMP", "Ecocert", "COSMOS", "REACH",
    "Halal", "Kosher", "Fair Trade", "USDA Organic", "RSB", "RSPO",
]
AVAILABILITY_STATES = {"in_stock", "limited", "out_of_stock"}


# ── Dataclasses ───────────────────────────────────────────────────────────────

@dataclass
class SupplierProfile:
    id:             str
    name:           str
    email:          str
    country:        str
    certifications: List[str]
    status:         str    # "pending" | "active" | "suspended"
    api_key_hash:   str    # SHA-256 hex; empty until approved
    api_key_prefix: str    # first 12 chars; empty until approved
    joined_at:      str
    description:    str = ""


@dataclass
class SupplierListing:
    id:              str
    supplier_id:     str
    supplier_name:   str
    ingredient_name: str
    price_per_kg:    float
    currency:        str
    moq_kg:          float
    lead_time_days:  int
    availability:    str
    certifications:  List[str]
    valid_until:     Optional[str]
    created_at:      str
    updated_at:      str


@dataclass
class PricingResult:
    ingredient:        str
    pct:               float
    best_price_per_kg: float
    best_supplier:     Optional[str]
    best_supplier_id:  Optional[str]
    lead_time_days:    Optional[int]
    moq_kg:            Optional[float]
    source:            str    # "supplier_listing" | "db_static"
    n_suppliers:       int
    price_spread_pct:  Optional[float]


@dataclass
class IngredientRisk:
    ingredient:       str
    pct:              float
    n_suppliers:      int
    single_sourced:   bool
    countries:        List[str]
    best_price:       Optional[float]
    price_spread_pct: Optional[float]
    min_lead_days:    Optional[int]
    max_lead_days:    Optional[int]
    risk_level:       str    # "low" | "medium" | "high"


@dataclass
class SupplyRiskReport:
    ingredient_risks:     List[IngredientRisk]
    overall_risk:         str
    single_sourced_count: int
    uncovered_count:      int
    geo_concentration:    Dict[str, float]   # country → % of blend weight
    weighted_lead_time:   Optional[float]
    total_ingredients:    int


# ── Persistence ───────────────────────────────────────────────────────────────

def _load_suppliers() -> List[SupplierProfile]:
    if not os.path.exists(_SUPPLIERS_PATH):
        return []
    try:
        with open(_SUPPLIERS_PATH) as f:
            return [SupplierProfile(**r) for r in json.load(f)]
    except Exception:
        return []


def _save_suppliers(suppliers: List[SupplierProfile]) -> None:
    try:
        with open(_SUPPLIERS_PATH, "w") as f:
            json.dump([asdict(s) for s in suppliers], f, indent=2)
    except Exception:
        pass


def _load_listings() -> List[SupplierListing]:
    if not os.path.exists(_LISTINGS_PATH):
        return []
    try:
        with open(_LISTINGS_PATH) as f:
            return [SupplierListing(**r) for r in json.load(f)]
    except Exception:
        return []


def _save_listings(listings: List[SupplierListing]) -> None:
    try:
        with open(_LISTINGS_PATH, "w") as f:
            json.dump([asdict(l) for l in listings], f, indent=2)
    except Exception:
        pass


# ── Supplier lifecycle ────────────────────────────────────────────────────────

def register_supplier(
    name: str,
    email: str,
    country: str,
    certifications: Optional[List[str]] = None,
    description: str = "",
) -> SupplierProfile:
    """Register a new supplier (status=pending). No key issued until approved."""
    suppliers = _load_suppliers()
    if any(s.email.lower() == email.lower() for s in suppliers):
        raise ValueError(f"A supplier with email '{email}' is already registered.")
    profile = SupplierProfile(
        id=secrets.token_hex(6),
        name=name,
        email=email,
        country=country,
        certifications=certifications or [],
        status="pending",
        api_key_hash="",
        api_key_prefix="",
        joined_at=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        description=description,
    )
    suppliers.append(profile)
    _save_suppliers(suppliers)
    return profile


def approve_supplier(supplier_id: str) -> Tuple[Optional[SupplierProfile], str]:
    """
    Approve a pending supplier and issue an ifs_ key.
    Returns (profile, raw_key). raw_key must be shown to the supplier exactly once.
    """
    suppliers = _load_suppliers()
    for s in suppliers:
        if s.id == supplier_id:
            raw_key = "ifs_" + secrets.token_urlsafe(32)
            s.api_key_hash   = hashlib.sha256(raw_key.encode()).hexdigest()
            s.api_key_prefix = raw_key[:12]
            s.status         = "active"
            _save_suppliers(suppliers)
            return s, raw_key
    return None, ""


def suspend_supplier(supplier_id: str) -> bool:
    suppliers = _load_suppliers()
    for s in suppliers:
        if s.id == supplier_id:
            s.status = "suspended"
            _save_suppliers(suppliers)
            return True
    return False


def validate_supplier_key(raw_key: str) -> Optional[SupplierProfile]:
    """Return the active SupplierProfile for a raw ifs_ key, or None."""
    if not raw_key:
        return None
    key_hash = hashlib.sha256(raw_key.encode()).hexdigest()
    for s in _load_suppliers():
        if s.api_key_hash == key_hash and s.status == "active":
            return s
    return None


def get_supplier(supplier_id: str) -> Optional[SupplierProfile]:
    return next((s for s in _load_suppliers() if s.id == supplier_id), None)


def list_suppliers(status: Optional[str] = None) -> List[SupplierProfile]:
    sups = _load_suppliers()
    return [s for s in sups if s.status == status] if status else sups


# ── Listing management ────────────────────────────────────────────────────────

def submit_listing(
    supplier_id:     str,
    ingredient_name: str,
    price_per_kg:    float,
    currency:        str = "USD",
    moq_kg:          float = 25.0,
    lead_time_days:  int = 14,
    availability:    str = "in_stock",
    certifications:  Optional[List[str]] = None,
    valid_until:     Optional[str] = None,
) -> SupplierListing:
    """Create or update a listing. Upserts on (supplier_id, ingredient_name)."""
    supplier = get_supplier(supplier_id)
    if not supplier:
        raise ValueError(f"Supplier '{supplier_id}' not found.")
    if supplier.status != "active":
        raise ValueError("Only active suppliers can submit listings.")
    if availability not in AVAILABILITY_STATES:
        raise ValueError(f"availability must be one of {AVAILABILITY_STATES}")
    if price_per_kg <= 0:
        raise ValueError("price_per_kg must be positive.")

    listings = _load_listings()
    now = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())

    for lst in listings:
        if lst.supplier_id == supplier_id and lst.ingredient_name.lower() == ingredient_name.lower():
            lst.price_per_kg   = price_per_kg
            lst.currency       = currency
            lst.moq_kg         = moq_kg
            lst.lead_time_days = lead_time_days
            lst.availability   = availability
            lst.certifications = certifications or []
            lst.valid_until    = valid_until
            lst.updated_at     = now
            _save_listings(listings)
            return lst

    listing = SupplierListing(
        id=secrets.token_hex(5),
        supplier_id=supplier_id,
        supplier_name=supplier.name,
        ingredient_name=ingredient_name,
        price_per_kg=price_per_kg,
        currency=currency,
        moq_kg=moq_kg,
        lead_time_days=lead_time_days,
        availability=availability,
        certifications=certifications or [],
        valid_until=valid_until,
        created_at=now,
        updated_at=now,
    )
    listings.append(listing)
    _save_listings(listings)
    return listing


def get_listings_for_supplier(supplier_id: str) -> List[SupplierListing]:
    return [l for l in _load_listings() if l.supplier_id == supplier_id]


def get_active_listings_for_ingredient(ingredient_name: str) -> List[SupplierListing]:
    name_lower = ingredient_name.lower()
    return [
        l for l in _load_listings()
        if l.ingredient_name.lower() == name_lower and l.availability != "out_of_stock"
    ]


def delete_listing(supplier_id: str, listing_id: str) -> bool:
    listings = _load_listings()
    filtered = [l for l in listings if not (l.id == listing_id and l.supplier_id == supplier_id)]
    if len(filtered) < len(listings):
        _save_listings(filtered)
        return True
    return False


# ── Pricing optimizer ─────────────────────────────────────────────────────────

def get_best_pricing(blend: Dict[str, float], db=None) -> List[PricingResult]:
    """
    For each ingredient in the blend return the lowest in-stock price.
    Falls back to db static Cost_USD_kg when no supplier listing exists.
    """
    results = []
    for ingredient, pct in blend.items():
        active = get_active_listings_for_ingredient(ingredient)

        if active:
            best = min(active, key=lambda l: l.price_per_kg)
            prices = [l.price_per_kg for l in active]
            spread = ((max(prices) - min(prices)) / min(prices) * 100) if len(prices) > 1 else None
            results.append(PricingResult(
                ingredient=ingredient,
                pct=pct,
                best_price_per_kg=best.price_per_kg,
                best_supplier=best.supplier_name,
                best_supplier_id=best.supplier_id,
                lead_time_days=best.lead_time_days,
                moq_kg=best.moq_kg,
                source="supplier_listing",
                n_suppliers=len({l.supplier_id for l in active}),
                price_spread_pct=round(spread, 1) if spread is not None else None,
            ))
        else:
            static_price = None
            if db is not None:
                try:
                    row = db[db["Ingredient"].str.lower() == ingredient.lower()]
                    if not row.empty:
                        static_price = float(row.iloc[0]["Cost_USD_kg"])
                except Exception:
                    pass
            results.append(PricingResult(
                ingredient=ingredient,
                pct=pct,
                best_price_per_kg=static_price or 0.0,
                best_supplier=None,
                best_supplier_id=None,
                lead_time_days=None,
                moq_kg=None,
                source="db_static",
                n_suppliers=0,
                price_spread_pct=None,
            ))
    return results


# ── Risk scoring ──────────────────────────────────────────────────────────────

def score_supply_risk(blend: Dict[str, float], db=None) -> SupplyRiskReport:
    """
    Risk levels per ingredient:
      high   — 0 suppliers (uncovered) or 1 supplier single-sourced
      medium — 2 suppliers, or all suppliers in one country
      low    — 3+ suppliers across 2+ countries
    """
    ingredient_risks: List[IngredientRisk] = []
    geo_weight: Dict[str, float] = {}
    total_pct = sum(blend.values()) or 1.0
    lead_weighted: List[float] = []

    for ingredient, pct in blend.items():
        active = get_active_listings_for_ingredient(ingredient)
        supplier_ids = list({l.supplier_id for l in active})
        countries = [
            get_supplier(sid).country
            for sid in supplier_ids
            if get_supplier(sid) is not None
        ]
        unique_countries = list(set(countries))
        n_sup = len(supplier_ids)
        prices = [l.price_per_kg for l in active]
        leads  = [l.lead_time_days for l in active]
        spread = ((max(prices) - min(prices)) / min(prices) * 100) if len(prices) > 1 else None

        if n_sup == 0:
            risk_level = "high"
        elif n_sup == 1 or len(unique_countries) == 1:
            risk_level = "medium"
        else:
            risk_level = "low"

        ingredient_risks.append(IngredientRisk(
            ingredient=ingredient,
            pct=pct,
            n_suppliers=n_sup,
            single_sourced=(n_sup == 1),
            countries=unique_countries,
            best_price=min(prices) if prices else None,
            price_spread_pct=round(spread, 1) if spread is not None else None,
            min_lead_days=min(leads) if leads else None,
            max_lead_days=max(leads) if leads else None,
            risk_level=risk_level,
        ))

        for c in unique_countries:
            geo_weight[c] = geo_weight.get(c, 0.0) + (pct / total_pct * 100)

        if leads:
            lead_weighted.append(min(leads) * (pct / total_pct))

    single_sourced_count = sum(1 for r in ingredient_risks if r.single_sourced)
    uncovered_count      = sum(1 for r in ingredient_risks if r.n_suppliers == 0)
    high_count           = sum(1 for r in ingredient_risks if r.risk_level == "high")
    medium_count         = sum(1 for r in ingredient_risks if r.risk_level == "medium")

    if high_count > 1 or uncovered_count > 1:
        overall_risk = "high"
    elif high_count == 1 or medium_count > 0:
        overall_risk = "medium"
    else:
        overall_risk = "low"

    return SupplyRiskReport(
        ingredient_risks=ingredient_risks,
        overall_risk=overall_risk,
        single_sourced_count=single_sourced_count,
        uncovered_count=uncovered_count,
        geo_concentration={k: round(v, 1) for k, v in geo_weight.items()},
        weighted_lead_time=round(sum(lead_weighted), 1) if lead_weighted else None,
        total_ingredients=len(blend),
    )


# ── Demand signals ────────────────────────────────────────────────────────────

def get_demand_signals(supplier_id: str) -> List[Dict]:
    """
    How many saved formulations use each ingredient this supplier lists.
    Reads from persistence (memory store + Supabase when available).
    """
    from modules.persistence import _MEMORY_STORE, load_all_projects
    listings = get_listings_for_supplier(supplier_id)
    if not listings:
        return []

    ingredient_names = {l.ingredient_name.lower() for l in listings}
    try:
        projects = load_all_projects(limit=500)
    except Exception:
        projects = _MEMORY_STORE.get("projects", [])

    counts: Dict[str, int] = {n: 0 for n in ingredient_names}
    for proj in projects:
        for ing in (proj.get("blend") or {}):
            key = ing.lower()
            if key in counts:
                counts[key] += 1

    signals = []
    for listing in listings:
        signals.append({
            "ingredient":        listing.ingredient_name,
            "listing_id":        listing.id,
            "formulation_count": counts.get(listing.ingredient_name.lower(), 0),
            "availability":      listing.availability,
            "price_per_kg":      listing.price_per_kg,
            "currency":          listing.currency,
        })
    signals.sort(key=lambda x: x["formulation_count"], reverse=True)
    return signals
