"""
modules/pubchem.py
PubChem chemical identity enrichment for IntelliForm.

Uses pubchempy (PubChem PUG REST API) to fetch CAS numbers, SMILES,
molecular weight, and molecular formula for ingredient names.

Free public API — no key required.
Reference: Kim et al., J Cheminformatics 2015, doi:10.1186/s13321-015-0060-4
"""
from __future__ import annotations

import re
import time
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from typing import Dict, Optional

logger = logging.getLogger(__name__)

_CACHE: Dict[str, tuple] = {}
_CACHE_TTL = 7 * 24 * 3600   # 7 days
_MAX_WORKERS = 8

_CAS_RE = re.compile(r'^\d{2,7}-\d{2}-\d$')


@dataclass
class PubChemResult:
    name: str
    cid: Optional[int] = None
    cas: Optional[str] = None
    smiles: Optional[str] = None
    iupac_name: Optional[str] = None
    molecular_formula: Optional[str] = None
    molecular_weight: Optional[float] = None
    inchikey: Optional[str] = None
    error: Optional[str] = None


def lookup(name: str) -> PubChemResult:
    """Look up one ingredient by name. Results are in-process cached for 7 days."""
    now = time.time()
    if name in _CACHE:
        ts, result = _CACHE[name]
        if now - ts < _CACHE_TTL:
            return result

    try:
        import pubchempy as pcp

        compounds = pcp.get_compounds(name, "name")
        if not compounds:
            result = PubChemResult(name=name, error="not found")
            _CACHE[name] = (now, result)
            return result

        c = compounds[0]

        cas = next((s for s in (c.synonyms or []) if _CAS_RE.match(s)), None)

        mw = None
        try:
            mw = float(c.molecular_weight) if c.molecular_weight else None
        except (TypeError, ValueError):
            pass

        result = PubChemResult(
            name=name,
            cid=c.cid,
            cas=cas,
            smiles=c.isomeric_smiles,
            iupac_name=c.iupac_name,
            molecular_formula=c.molecular_formula,
            molecular_weight=mw,
            inchikey=c.inchikey,
        )

    except Exception as exc:
        logger.debug(f"[pubchem] lookup failed for {name!r}: {exc}")
        result = PubChemResult(name=name, error=str(exc))

    _CACHE[name] = (now, result)
    return result


def enrich_blend(blend: dict) -> Dict[str, dict]:
    """
    Parallel PubChem lookup for all ingredients in a blend.

    Args:
        blend: {ingredient_name: weight_percent}

    Returns:
        {ingredient_name: PubChemResult fields as dict}
    """
    names = list(blend.keys())
    results: Dict[str, PubChemResult] = {}

    with ThreadPoolExecutor(max_workers=min(_MAX_WORKERS, len(names) or 1)) as ex:
        futures = {ex.submit(lookup, name): name for name in names}
        for future in as_completed(futures):
            name = futures[future]
            try:
                results[name] = future.result()
            except Exception as exc:
                results[name] = PubChemResult(name=name, error=str(exc))

    return {n: r.__dict__ for n, r in results.items()}
