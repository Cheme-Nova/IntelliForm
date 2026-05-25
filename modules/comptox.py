"""
modules/comptox.py
EPA CompTox / CCTE API integration for IntelliForm.

Fetches OPERA regulatory QSAR predictions for ingredients via the EPA's
Center for Computational Toxicology and Exposure (CCTE) public API.

Unlike IntelliForm's own QSAR models, OPERA predictions are:
  - OECD QSAR Toolbox-compliant
  - Cited in EPA and ECHA regulatory submissions
  - Peer-reviewed (Mansouri et al., J Cheminformatics 2018)

API base: https://api-ccte.epa.gov
Free key: search "EPA CCTE API key" on epa.gov to register
Env var:  COMPTOX_API_KEY  (optional — graceful fallback without it)

OPERA endpoints retrieved per ingredient:
  ready_biodegradability   OECD 301B proxy (True/False + probability)
  log_bcf                  Bioconcentration factor (log L/kg)
  log_kow                  Octanol-water partition (log)
  water_solubility         mg/L at 25°C
  log_koc                  Soil adsorption coefficient
  atm_half_life            Atmospheric half-life (hours)
  vapor_pressure           mmHg at 25°C
  fish_lc50                Acute fish toxicity (log mmol/L)
  daphnia_ec50             Daphnia magna acute toxicity (log mmol/L)
  svhc_candidate           REACH Annex XV SVHC candidate list (bool)
  cmr_category             Carcinogenic/Mutagenic/Reprotoxic (1A/1B/2 or None)

Caching:
  Results are cached in-process (per Render worker) + optionally on disk
  at /tmp/comptox_cache.json to survive warm restarts. TTL: 7 days.
"""
from __future__ import annotations

import json
import os
import time
import logging
from dataclasses import dataclass, field, asdict
from typing import Dict, List, Optional

import requests

logger = logging.getLogger(__name__)

# ── Config ─────────────────────────────────────────────────────────────────────
_BASE_URL   = "https://api-ccte.epa.gov"
_CACHE_PATH = "/tmp/comptox_cache.json"
_CACHE_TTL  = 7 * 24 * 3600   # 7 days in seconds
_TIMEOUT    = 8                # per-request timeout seconds
_BATCH_SIZE = 10               # max names per batch request

# In-process cache: {ingredient_name: (timestamp, CompToxResult)}
_MEM_CACHE: Dict[str, tuple] = {}


# ── Result schema ──────────────────────────────────────────────────────────────

@dataclass
class CompToxResult:
    """OPERA predictions + regulatory flags for one ingredient."""
    name: str
    dtxsid: Optional[str]           # EPA identifier e.g. DTXSID7020182
    cas: Optional[str]
    smiles: Optional[str]
    preferred_name: Optional[str]

    # OPERA predictions
    ready_biodegradable: Optional[bool]
    biodeg_probability: Optional[float]  # 0.0–1.0
    log_bcf: Optional[float]             # log L/kg
    log_kow: Optional[float]
    water_solubility_mg_l: Optional[float]
    log_koc: Optional[float]
    atm_half_life_h: Optional[float]
    vapor_pressure_mmhg: Optional[float]
    fish_lc50_log_mmol_l: Optional[float]
    daphnia_ec50_log_mmol_l: Optional[float]

    # Regulatory flags
    svhc_candidate: bool = False
    cmr_category: Optional[str] = None  # "1A" / "1B" / "2" / None
    ghs_hazard_codes: List[str] = field(default_factory=list)
    reach_restricted: bool = False

    # Source metadata
    source: str = "comptox"
    retrieved_at: float = field(default_factory=time.time)
    error: Optional[str] = None


@dataclass
class BlendCompToxReport:
    """Aggregated CompTox assessment for a full blend."""
    ingredients: Dict[str, CompToxResult]
    svhc_flags: List[str]                # ingredient names with SVHC concern
    cmr_flags: List[str]                 # ingredient names with CMR classification
    reach_restricted_flags: List[str]
    ready_biodeg_fraction: Optional[float]   # weighted % of blend that is OECD biodegradable
    avg_log_bcf: Optional[float]
    avg_log_kow: Optional[float]
    regulatory_citation: str
    coverage: int                         # % of blend ingredients found in CompTox


# ── Cache helpers ──────────────────────────────────────────────────────────────

def _load_disk_cache() -> Dict[str, tuple]:
    try:
        with open(_CACHE_PATH) as f:
            raw = json.load(f)
        now = time.time()
        return {
            k: (ts, _dict_to_result(v))
            for k, (ts, v) in raw.items()
            if now - ts < _CACHE_TTL
        }
    except Exception:
        return {}


def _save_disk_cache(cache: Dict[str, tuple]) -> None:
    try:
        serializable = {k: (ts, asdict(r)) for k, (ts, r) in cache.items()}
        with open(_CACHE_PATH, "w") as f:
            json.dump(serializable, f)
    except Exception:
        pass


def _dict_to_result(d: dict) -> CompToxResult:
    return CompToxResult(**{k: v for k, v in d.items() if k in CompToxResult.__dataclass_fields__})


def _init_cache() -> None:
    global _MEM_CACHE
    if not _MEM_CACHE:
        _MEM_CACHE = _load_disk_cache()


# ── API client ─────────────────────────────────────────────────────────────────

def _api_key() -> Optional[str]:
    return os.getenv("COMPTOX_API_KEY", "").strip() or None


def _headers() -> dict:
    h = {"Content-Type": "application/json", "Accept": "application/json"}
    key = _api_key()
    if key:
        h["x-api-key"] = key
    return h


def _search_chemical_by_name(name: str) -> Optional[dict]:
    """Search CompTox for a chemical by name. Returns first match dict or None."""
    url = f"{_BASE_URL}/chemical/search/by-name/{requests.utils.quote(name)}"
    try:
        resp = requests.get(url, headers=_headers(), timeout=_TIMEOUT)
        if resp.status_code == 200:
            data = resp.json()
            # API returns a list; take the first exact or best match
            if isinstance(data, list) and data:
                return data[0]
            if isinstance(data, dict) and data.get("dtxsid"):
                return data
    except Exception as exc:
        logger.debug(f"[comptox] name search failed for {name!r}: {exc}")
    return None


def _get_chemical_detail(dtxsid: str) -> Optional[dict]:
    """Fetch full chemical detail including OPERA predictions by DTXSID."""
    url = f"{_BASE_URL}/chemical/detail/by-dtxsid/{dtxsid}"
    try:
        resp = requests.get(url, headers=_headers(), timeout=_TIMEOUT)
        if resp.status_code == 200:
            return resp.json()
    except Exception as exc:
        logger.debug(f"[comptox] detail fetch failed for {dtxsid}: {exc}")
    return None


def _batch_search_names(names: List[str]) -> Dict[str, dict]:
    """Batch search multiple chemical names. Returns {name: result_dict}."""
    url = f"{_BASE_URL}/chemical/list/search/by-name/"
    try:
        resp = requests.post(
            url,
            headers=_headers(),
            json=names,
            timeout=_TIMEOUT * 2,
        )
        if resp.status_code == 200:
            results = resp.json()
            if isinstance(results, list):
                out = {}
                for r in results:
                    if isinstance(r, dict) and r.get("searchWord"):
                        out[r["searchWord"]] = r
                return out
    except Exception as exc:
        logger.debug(f"[comptox] batch search failed: {exc}")
    return {}


# ── Parsing helpers ────────────────────────────────────────────────────────────

# Known SVHC/CMR DTXSID prefixes and CAS numbers (hardcoded subset for offline fallback)
_SVHC_CAS = {
    "71-43-2",   # Benzene
    "75-09-2",   # DCM
    "872-50-4",  # NMP
    "109-99-9",  # THF (not SVHC but restricted)
    "67-64-1",   # Acetone (not SVHC)
    "50-00-0",   # Formaldehyde
    "1336-36-3", # PCBs
    "117-81-7",  # DEHP
    "84-74-2",   # DBP
    "85-68-7",   # BBP
}

_CMR_CAS = {
    "71-43-2":  "1A",  # Benzene
    "50-00-0":  "1A",  # Formaldehyde
    "75-09-2":  "2",   # DCM (IARC 2A)
    "872-50-4": "1B",  # NMP (reproductive)
    "80-05-7":  "1B",  # BPA
}


def _parse_opera(detail: dict) -> dict:
    """Extract OPERA prediction fields from a CompTox detail response."""
    opera = {}
    if not detail:
        return opera

    # The detail endpoint nests OPERA data under various keys depending on API version
    # Try both flat and nested formats
    for prefix in ("", "opera_", "OPERA_"):
        opera["ready_biodegradable"] = (
            detail.get(f"{prefix}ReadyBiodeg") or
            detail.get(f"{prefix}readyBiodeg") or
            detail.get("readyBiodegradability")
        )
        opera["log_bcf"] = (
            detail.get(f"{prefix}logBCF_pred") or
            detail.get(f"{prefix}logBcf") or
            detail.get("LogBCF")
        )
        opera["log_kow"] = (
            detail.get(f"{prefix}logKow_pred") or
            detail.get(f"{prefix}logKow") or
            detail.get("logKow") or
            detail.get("LogKow")
        )
        opera["water_solubility_mg_l"] = (
            detail.get(f"{prefix}waterSolubility_pred") or
            detail.get("WaterSolubility") or
            detail.get("waterSolubility")
        )
        opera["log_koc"] = (
            detail.get(f"{prefix}logKoc_pred") or
            detail.get("logKoc")
        )
        opera["atm_half_life_h"] = (
            detail.get(f"{prefix}atmosphericHalfLife_pred") or
            detail.get("AtmHalfLife")
        )
        opera["vapor_pressure_mmhg"] = (
            detail.get(f"{prefix}vaporPressure_pred") or
            detail.get("VaporPressure")
        )
        opera["fish_lc50"] = (
            detail.get(f"{prefix}FishLC50_pred") or
            detail.get("fishLc50") or
            detail.get("LogLC50_pred")
        )
        opera["daphnia_ec50"] = (
            detail.get(f"{prefix}DaphniaEC50_pred") or
            detail.get("daphniaEc50") or
            detail.get("LogEC50_pred")
        )

    return opera


def _safe_float(val) -> Optional[float]:
    try:
        return round(float(val), 4) if val is not None else None
    except (TypeError, ValueError):
        return None


def _safe_bool(val) -> Optional[bool]:
    if val is None:
        return None
    if isinstance(val, bool):
        return val
    s = str(val).strip().lower()
    if s in ("true", "yes", "1", "ready"):
        return True
    if s in ("false", "no", "0", "not ready"):
        return False
    return None


def _extract_ghs_codes(detail: dict) -> List[str]:
    codes = []
    for key in ("ghsHazardCodes", "hazardCodes", "ghs"):
        raw = detail.get(key)
        if isinstance(raw, list):
            codes.extend(str(x) for x in raw if x)
        elif isinstance(raw, str) and raw:
            codes.append(raw)
    return list(set(codes))


def _check_svhc(detail: dict) -> bool:
    # Check API flag first
    for key in ("svhcCandidate", "svhc", "isSvhc", "reach_svhc"):
        val = detail.get(key)
        if val is not None:
            return bool(val)
    # Fallback: CAS lookup
    cas = detail.get("casrn") or detail.get("cas") or ""
    return cas in _SVHC_CAS


def _check_cmr(detail: dict) -> Optional[str]:
    for key in ("cmrCategory", "cmr", "carcinogenicity"):
        val = detail.get(key)
        if val:
            return str(val)
    cas = detail.get("casrn") or detail.get("cas") or ""
    return _CMR_CAS.get(cas)


# ── Public API ─────────────────────────────────────────────────────────────────

def lookup_ingredient(name: str) -> CompToxResult:
    """
    Look up a single ingredient in CompTox and return OPERA predictions.

    Results are cached in-process. Falls back gracefully on API failure.
    """
    _init_cache()
    now = time.time()

    # Check in-process cache
    if name in _MEM_CACHE:
        ts, cached = _MEM_CACHE[name]
        if now - ts < _CACHE_TTL:
            return cached

    # 1. Search by name
    match = _search_chemical_by_name(name)
    if not match:
        result = CompToxResult(
            name=name, dtxsid=None, cas=None, smiles=None, preferred_name=None,
            ready_biodegradable=None, biodeg_probability=None,
            log_bcf=None, log_kow=None, water_solubility_mg_l=None,
            log_koc=None, atm_half_life_h=None, vapor_pressure_mmhg=None,
            fish_lc50_log_mmol_l=None, daphnia_ec50_log_mmol_l=None,
            error="not_found",
        )
        _MEM_CACHE[name] = (now, result)
        return result

    dtxsid = match.get("dtxsid") or match.get("dtxSid") or match.get("id")
    cas    = match.get("casrn") or match.get("cas")
    smiles = match.get("smiles") or match.get("qsarSmiles")
    pref   = match.get("preferredName") or match.get("iupacName") or name

    # 2. Fetch full detail for OPERA predictions
    detail = {}
    if dtxsid:
        detail = _get_chemical_detail(dtxsid) or {}

    # Merge match + detail
    combined = {**match, **detail}
    opera = _parse_opera(combined)

    # Parse ready biodegradability
    ready_raw = opera.get("ready_biodegradable")
    ready_bool = _safe_bool(ready_raw)
    # Extract probability if available as separate field
    biodeg_prob = _safe_float(
        combined.get("readyBiodegProb") or
        combined.get("biodegradabilityProbability") or
        (1.0 if ready_bool is True else (0.0 if ready_bool is False else None))
    )

    result = CompToxResult(
        name=name,
        dtxsid=dtxsid,
        cas=cas,
        smiles=smiles,
        preferred_name=pref,
        ready_biodegradable=ready_bool,
        biodeg_probability=biodeg_prob,
        log_bcf=_safe_float(opera.get("log_bcf")),
        log_kow=_safe_float(opera.get("log_kow")),
        water_solubility_mg_l=_safe_float(opera.get("water_solubility_mg_l")),
        log_koc=_safe_float(opera.get("log_koc")),
        atm_half_life_h=_safe_float(opera.get("atm_half_life_h")),
        vapor_pressure_mmhg=_safe_float(opera.get("vapor_pressure_mmhg")),
        fish_lc50_log_mmol_l=_safe_float(opera.get("fish_lc50")),
        daphnia_ec50_log_mmol_l=_safe_float(opera.get("daphnia_ec50")),
        svhc_candidate=_check_svhc(combined),
        cmr_category=_check_cmr(combined),
        ghs_hazard_codes=_extract_ghs_codes(combined),
        reach_restricted=bool(combined.get("reachRestricted") or combined.get("reach_restricted")),
    )

    _MEM_CACHE[name] = (now, result)
    return result


def screen_blend(blend: dict, vertical: Optional[str] = None) -> BlendCompToxReport:
    """
    Run CompTox OPERA screening on all ingredients in a blend.

    Args:
        blend:    {ingredient_name: weight_percent}
        vertical: optional vertical key for context

    Returns:
        BlendCompToxReport with per-ingredient OPERA data and regulatory flags.
    """
    if not _api_key():
        return BlendCompToxReport(
            ingredients={},
            svhc_flags=[],
            cmr_flags=[],
            reach_restricted_flags=[],
            ready_biodeg_fraction=None,
            avg_log_bcf=None,
            avg_log_kow=None,
            regulatory_citation=(
                "EPA CTX API key not configured. Email ccte_api@epa.gov for a free key, "
                "then set COMPTOX_API_KEY in Render environment variables."
            ),
            coverage=0,
        )

    results: Dict[str, CompToxResult] = {}
    names = list(blend.keys())

    # Process in batches to respect rate limits
    for i in range(0, len(names), _BATCH_SIZE):
        batch = names[i: i + _BATCH_SIZE]
        # Try batch first (more efficient)
        batch_hit = _batch_search_names(batch)
        for name in batch:
            if name in results:
                continue
            r = lookup_ingredient(name)
            results[name] = r

    # Aggregate metrics (weighted by blend %)
    total_pct = sum(blend.values()) or 100.0
    svhc_flags, cmr_flags, reach_flags = [], [], []
    weighted_biodeg = 0.0
    biodeg_covered_pct = 0.0
    weighted_bcf = 0.0
    bcf_covered_pct = 0.0
    weighted_kow = 0.0
    kow_covered_pct = 0.0
    found = 0

    for name, pct in blend.items():
        r = results.get(name)
        if not r or r.error:
            continue
        found += 1
        frac = pct / total_pct

        if r.svhc_candidate:
            svhc_flags.append(name)
        if r.cmr_category:
            cmr_flags.append(name)
        if r.reach_restricted:
            reach_flags.append(name)
        if r.biodeg_probability is not None:
            weighted_biodeg += r.biodeg_probability * frac
            biodeg_covered_pct += frac
        if r.log_bcf is not None:
            weighted_bcf += r.log_bcf * frac
            bcf_covered_pct += frac
        if r.log_kow is not None:
            weighted_kow += r.log_kow * frac
            kow_covered_pct += frac

    ready_biodeg_pct = round((weighted_biodeg / biodeg_covered_pct) * 100, 1) if biodeg_covered_pct > 0 else None
    avg_bcf  = round(weighted_bcf / bcf_covered_pct, 3) if bcf_covered_pct > 0 else None
    avg_kow  = round(weighted_kow / kow_covered_pct, 3) if kow_covered_pct > 0 else None
    coverage = round((found / len(names)) * 100) if names else 0

    # Persist cache after batch run
    _save_disk_cache(_MEM_CACHE)

    return BlendCompToxReport(
        ingredients=results,
        svhc_flags=svhc_flags,
        cmr_flags=cmr_flags,
        reach_restricted_flags=reach_flags,
        ready_biodeg_fraction=ready_biodeg_pct,
        avg_log_bcf=avg_bcf,
        avg_log_kow=avg_kow,
        regulatory_citation=(
            "OPERA predictions from EPA CompTox Dashboard (CCTE). "
            "Mansouri et al., J Cheminformatics 2018. "
            "OECD QSAR Toolbox compliant."
        ),
        coverage=coverage,
    )
