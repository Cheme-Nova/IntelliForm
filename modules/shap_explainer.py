"""
modules/shap_explainer.py
SHAP ingredient attribution for IntelliForm blend formulations.

For each ingredient in a blend:
  - TreeExplainer SHAP values on the active GBR model (Biodegradability, Ecotoxicity, Performance)
  - Weighted blend contribution vs blend average (signed, in endpoint units)
  - Top-N named chemical features driving each endpoint prediction

Falls back to weighted-delta attribution when SHAP is not installed (pip install shap).
TreeExplainer operates on the GBR estimator inside the sklearn Pipeline, using the
imputer+scaler-preprocessed feature matrix so indices match the original feature names.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

try:
    import shap as _shap
    SHAP_OK = True
except ImportError:
    SHAP_OK = False

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    _RDKIT_AVAIL = True
except ImportError:
    _RDKIT_AVAIL = False

from modules.qsar import (
    get_gbr_models,
    _smiles_to_features,
    _feature_names,
    predict_properties,
    MORGAN_NBITS,
)

# Physicochemical descriptor names appended after Morgan FP bits
_PHYS_DESCS = ["MW", "LogP", "TPSA", "HBondAcc", "HBondDon", "RotBonds", "FractionCSP3"]
_TARGETS = ["Biodegradability", "Ecotoxicity", "Performance"]


# ── Dataclasses ───────────────────────────────────────────────────────────────

@dataclass
class IngredientAttribution:
    ingredient: str
    pct: float
    bio_prediction: float
    etox_prediction: float
    perf_prediction: float
    # Pct-weighted delta from blend average (signed, in endpoint units)
    bio_contribution: float
    etox_contribution: float
    perf_contribution: float
    # Top N (feature_name, shap_value) pairs per endpoint
    bio_top_features:  List[Tuple[str, float]] = field(default_factory=list)
    etox_top_features: List[Tuple[str, float]] = field(default_factory=list)
    perf_top_features: List[Tuple[str, float]] = field(default_factory=list)
    shap_available: bool = False


@dataclass
class BlendAttribution:
    ingredients: List[IngredientAttribution]
    blend_bio:  float
    blend_etox: float
    blend_perf: float
    shap_available: bool
    method: str   # "shap_treexplainer" | "weighted_delta"


# ── Feature-naming helpers ────────────────────────────────────────────────────

def _bit_to_name(mol, bit_idx: int, bit_info: dict) -> str:
    """Decode a Morgan FP bit index to a human-readable fragment SMILES or descriptor name."""
    if bit_idx >= MORGAN_NBITS:
        offset = bit_idx - MORGAN_NBITS
        return _PHYS_DESCS[offset] if offset < len(_PHYS_DESCS) else f"desc_{offset}"
    if not _RDKIT_AVAIL or mol is None or bit_idx not in bit_info:
        return f"bit_{bit_idx}"
    try:
        atom_idx, radius = bit_info[bit_idx][0]
        if radius == 0:
            return f"{mol.GetAtomWithIdx(atom_idx).GetSymbol()}-center"
        env = Chem.FindAtomEnvironmentOfRadiusN(mol, radius, atom_idx)
        if not env:
            return f"bit_{bit_idx}"
        atom_set: set = set()
        for bond_idx in env:
            b = mol.GetBondWithIdx(bond_idx)
            atom_set.add(b.GetBeginAtomIdx())
            atom_set.add(b.GetEndAtomIdx())
        frag = AllChem.MolFragmentToSmiles(
            mol,
            atomsToUse=list(atom_set),
            bondsToUse=list(env),
            isomericSmiles=False,
        )
        if frag:
            return frag[:24]
        return f"env_r{radius}"
    except Exception:
        return f"bit_{bit_idx}"


def _top_features(
    shap_vals: np.ndarray,
    feat_names: List[str],
    mol,
    bit_info: dict,
    n: int,
    used_mordred: bool,
) -> List[Tuple[str, float]]:
    """Return top-n features by |SHAP value| with human-readable names."""
    indices = np.argsort(np.abs(shap_vals))[::-1][:n]
    out = []
    for idx in indices:
        val = float(shap_vals[idx])
        if abs(val) < 1e-7:
            continue
        if used_mordred:
            name = feat_names[idx] if idx < len(feat_names) else f"feat_{idx}"
        else:
            name = _bit_to_name(mol, int(idx), bit_info)
        out.append((name, round(val, 4)))
    return out


# ── SHAP computation ──────────────────────────────────────────────────────────

def _shap_for_pipe(pipe, X_2d: np.ndarray) -> Optional[np.ndarray]:
    """
    Run SHAP TreeExplainer on the GBR stage of a sklearn Pipeline.
    Preprocessing (imputer + scaler) is applied first so feature indices
    remain aligned with the original feature name list.
    Returns SHAP values of shape (n_features,), or None on any error.
    """
    if not SHAP_OK:
        return None
    try:
        X_pre = pipe[:-1].transform(X_2d)
        gbr   = pipe.named_steps["model"]
        explainer = _shap.TreeExplainer(gbr, feature_perturbation="tree_path_dependent")
        sv = explainer.shap_values(X_pre)
        # GBR regressor: sv is (n_samples, n_features) for 2-D input
        return sv[0] if (isinstance(sv, np.ndarray) and sv.ndim == 2) else np.array(sv)
    except Exception:
        return None


# ── Public API ────────────────────────────────────────────────────────────────

def explain_blend(
    blend: Dict[str, float],
    db: pd.DataFrame,
    top_n: int = 5,
) -> BlendAttribution:
    """
    Compute per-ingredient SHAP attribution across all three endpoints
    for every ingredient in `blend`.

    Returns a BlendAttribution with:
      - per-ingredient bio/etox/perf contribution (pct-weighted delta from blend mean)
      - top-n SHAP feature names per endpoint per ingredient (when SHAP available)
      - blend-level weighted averages

    Falls back to weighted-delta attribution when SHAP is not installed.
    """
    gbr_models, used_mordred = get_gbr_models()
    feat_names  = _feature_names(used_mordred)
    shap_avail  = SHAP_OK and (gbr_models is not None)
    method      = "shap_treexplainer" if shap_avail else "weighted_delta"

    # ── Predict each ingredient ───────────────────────────────────────────────
    preds: Dict[str, Tuple[float, float, float]] = {}
    for ing in blend:
        rows   = db[db["Ingredient"] == ing] if db is not None else pd.DataFrame()
        smiles = str(rows.iloc[0].get("SMILES", "")) if not rows.empty else ""
        qp     = predict_properties(smiles if smiles and smiles != "nan" else "")
        preds[ing] = (qp.biodegradability, qp.ecotoxicity, qp.performance)

    total_pct  = sum(blend.values()) or 1.0
    blend_bio  = sum(blend[i] / total_pct * preds[i][0] for i in blend)
    blend_etox = sum(blend[i] / total_pct * preds[i][1] for i in blend)
    blend_perf = sum(blend[i] / total_pct * preds[i][2] for i in blend)

    # ── Per-ingredient attribution ────────────────────────────────────────────
    attributions: List[IngredientAttribution] = []

    for ing, pct in blend.items():
        bio, etox, perf = preds[ing]
        w = pct / total_pct

        bio_contrib  = round(w * (bio  - blend_bio),  2)
        etox_contrib = round(w * (etox - blend_etox), 3)
        perf_contrib = round(w * (perf - blend_perf), 2)

        bio_top = etox_top = perf_top = []

        if shap_avail:
            rows   = db[db["Ingredient"] == ing] if db is not None else pd.DataFrame()
            smiles = str(rows.iloc[0].get("SMILES", "")) if not rows.empty else ""

            mol: Optional[object] = None
            bit_info: dict = {}

            if smiles and smiles != "nan" and _RDKIT_AVAIL:
                mol = Chem.MolFromSmiles(smiles)
                if mol is not None and not used_mordred:
                    AllChem.GetMorganFingerprintAsBitVect(
                        mol, 2, nBits=MORGAN_NBITS, bitInfo=bit_info
                    )

            feat, _ = _smiles_to_features(smiles) if smiles and smiles != "nan" else (None, False)
            if feat is not None:
                X_2d = np.nan_to_num(feat, nan=0.0).reshape(1, -1)
                for target in _TARGETS:
                    pipe = gbr_models.get(target)
                    if pipe is None:
                        continue
                    sv = _shap_for_pipe(pipe, X_2d)
                    if sv is None:
                        continue
                    tops = _top_features(sv, feat_names, mol, bit_info, top_n, used_mordred)
                    if target == "Biodegradability":
                        bio_top = tops
                    elif target == "Ecotoxicity":
                        etox_top = tops
                    else:
                        perf_top = tops

        attributions.append(IngredientAttribution(
            ingredient=ing,
            pct=pct,
            bio_prediction=bio,
            etox_prediction=etox,
            perf_prediction=perf,
            bio_contribution=bio_contrib,
            etox_contribution=etox_contrib,
            perf_contribution=perf_contrib,
            bio_top_features=bio_top,
            etox_top_features=etox_top,
            perf_top_features=perf_top,
            shap_available=shap_avail,
        ))

    return BlendAttribution(
        ingredients=attributions,
        blend_bio=round(blend_bio, 1),
        blend_etox=round(blend_etox, 2),
        blend_perf=round(blend_perf, 1),
        shap_available=shap_avail,
        method=method,
    )
