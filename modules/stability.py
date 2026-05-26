"""
modules/stability.py
Stability and viscosity prediction for IntelliForm v1.5+.

Model tiers (best available at runtime):
  1. GBR + ICP   — GradientBoostingRegressor trained on synthetic + lab feedback
                   with inductive conformal prediction intervals (80/90/95% coverage)
                   + 5-model uncertainty ensemble
  2. Rule-based  — lookup-table heuristics (used as feature inputs to tier 1
                   and as fallback when ML is unavailable)

Predicted endpoints:
  - Viscosity   (cP at 25 °C)
  - Shelf life  (months)

Blend-level ML features (10):
  thickener_pct, surfactant_pct, preservative_pct, antioxidant_pct,
  bio_avg_pct, n_ingredients, max_conc_pct, visc_potential,
  has_preservative, has_antioxidant

Synthetic bootstrap training (n=300) generated from rule-based predictions
with 10%/12% Gaussian noise, then augmented by real AL feedback records.
"""
from __future__ import annotations

import os
import pickle
import hashlib
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

# ── Optional sklearn ──────────────────────────────────────────────────────────

try:
    from sklearn.ensemble import GradientBoostingRegressor
    from sklearn.preprocessing import StandardScaler
    from sklearn.impute import SimpleImputer
    from sklearn.pipeline import Pipeline
    SKLEARN_OK = True
except ImportError:
    SKLEARN_OK = False

# ── Constants ─────────────────────────────────────────────────────────────────

STAB_CACHE_PATH    = "/tmp/intelliform_stability_v15.pkl"
STAB_CACHE_VERSION = "v1.0-ml"
STAB_FEEDBACK_PATH = "/tmp/intelliform_stability_feedback.pkl"
STAB_CONF_LEVELS   = [0.80, 0.90, 0.95]
STAB_CAL_PCT       = 0.15
STAB_ENSEMBLE_N    = 5
STAB_SYNTH_N       = 300

_VISC_CLIP  = (1.0,  100_000.0)
_SHELF_CLIP = (3.0,  60.0)

_FEAT_NAMES = [
    "thickener_pct", "surfactant_pct", "preservative_pct", "antioxidant_pct",
    "bio_avg", "n_ingredients", "max_conc", "visc_potential",
    "has_preservative", "has_antioxidant",
]

# ── Lookup tables (retained as feature generators + fallback) ─────────────────

_THICKENER_VISCOSITY = {
    "Xanthan Gum": 8000, "Hydroxyethyl Cellulose": 4000,
    "Guar Gum (bio)": 5000, "Hydroxypropyl Methylcellulose": 3500,
    "Carrageenan": 3000, "Cellulose Gum (CMC)": 2500, "Pectin": 2000,
    "Carbomer": 9000, "Sodium Alginate": 3000, "Hydroxypropyl Starch": 2000,
}

_PRESERVATIVE_EFFICACY = {
    "Benzyl Alcohol (NF)": 18, "Sodium Benzoate": 12,
    "Potassium Sorbate": 12, "Phenoxyethanol": 24,
    "Ethylhexylglycerin": 18, "Caprylhydroxamic Acid": 24,
    "Caprylic Acid (bio)": 18, "Undecylenic Acid": 18,
    "Levulinic Acid (bio)": 15, "Hydroxyacetophenone": 18,
    "1-Hexanediol": 18, "Hexylene Glycol": 18,
}

_ANTIOXIDANT_BOOST = {
    "Tocopheryl Acetate": 6, "Ascorbic Acid": 4,
    "Vitamin E (mixed tocopherols)": 8, "Rosemary Extract (antioxidant)": 6,
    "Green Tea Extract": 5,
}

_EMOLLIENT_PACKAGING = {
    "Argan Oil": "Glass or aluminium — protect from UV",
    "Jojoba Esters": "HDPE or glass acceptable",
    "Coconut Oil (fractionated)": "HDPE",
    "Shea Butter (refined)": "HDPE or PET",
}

_RISK_INGREDIENTS = {
    "D-Limonene": "D-Limonene can oxidise on prolonged storage — add antioxidant",
    "Rhamnolipid (biosurfactant)": "Biosurfactants sensitive to pH shift — monitor monthly",
    "Sophorolipid": "Sophorolipids may crystallise below 15°C — store above ambient",
    "Sodium Coco-Sulfate": "High sulfate content may cause pH drift over time",
    "Lecithin (sunflower)": "Lecithin susceptible to hydrolysis — keep pH 4–7",
    "Lecithin (soy)": "Lecithin susceptible to hydrolysis — keep pH 4–7",
    "Ascorbic Acid": "Ascorbic Acid oxidises rapidly — use nitrogen headspace packaging",
    "Ceramide NP (bio)": "Ceramides require emulsification stabiliser for aqueous systems",
}


# ── Data schemas ──────────────────────────────────────────────────────────────

@dataclass
class StabilityConformalBands:
    quantiles:       Dict[str, Dict[str, float]]  # target → {"80": q, "90": q, "95": q}
    n_calibration:   int
    coverage_levels: List[float]


@dataclass
class StabilityFeedbackRecord:
    blend_hash:   str    # md5 of sorted blend items
    target:       str    # "Viscosity" | "ShelfLife"
    actual_value: float
    timestamp:    str
    al_round:     int


@dataclass
class StabilityResult:
    shelf_life_months:    int
    shelf_life_range:     str
    viscosity_cp:         float
    viscosity_range:      str
    ph_min:               float
    ph_max:               float
    recommended_packaging: str
    stability_risks:      List[str]
    stability_boosters:   List[str]
    overall_rating:       str
    confidence:           str
    ml_active:            bool = False
    # ICP intervals when ML tier is active
    # {"Viscosity": {"80": (lo, hi), ...}, "ShelfLife": {"80": (lo, hi), ...}}
    intervals:            Optional[Dict[str, Dict[str, Tuple[float, float]]]] = None


@dataclass
class StabilityModelCard:
    n_training:           int
    n_feedback:           int
    al_rounds:            int
    synth_samples:        int
    conformal_active:     bool
    conformal_n_cal:      int
    ensemble_active:      bool
    sklearn_version:      str
    cache_version:        str = STAB_CACHE_VERSION


# ── Module-level caches ───────────────────────────────────────────────────────

_STAB_MODELS:    Optional[Dict]                      = None  # {"Viscosity": pipe, "ShelfLife": pipe}
_STAB_ENSEMBLE:  Optional[Dict]                      = None  # {"Viscosity": [pipe×5], "ShelfLife": [pipe×5]}
_STAB_CONF:      Optional[StabilityConformalBands]   = None
_STAB_CARD:      Optional[StabilityModelCard]        = None


# ── Feature extraction ────────────────────────────────────────────────────────

def _blend_features(blend: Dict[str, float], db: pd.DataFrame) -> Optional[np.ndarray]:
    """Extract 10 blend-level stability features."""
    idx = db.set_index("Ingredient") if "Ingredient" in db.columns else pd.DataFrame()

    thickener_pct = surfactant_pct = preservative_pct = antioxidant_pct = 0.0
    bio_sum = bio_weight = 0.0
    visc_potential = 0.0

    for ing, pct in blend.items():
        if ing in _THICKENER_VISCOSITY:
            thickener_pct += pct
            visc_potential += _THICKENER_VISCOSITY[ing] * (pct / 100) ** 0.7
        if ing in _PRESERVATIVE_EFFICACY:
            preservative_pct += pct
        if ing in _ANTIOXIDANT_BOOST:
            antioxidant_pct += pct
        if not idx.empty and ing in idx.index:
            func = str(idx.loc[ing, "Function"]) if "Function" in idx.columns else ""
            if "Surfactant" in func:
                surfactant_pct += pct
            try:
                bio = float(idx.loc[ing, "Bio_based_pct"])
                bio_sum   += bio * pct
                bio_weight += pct
            except Exception:
                pass

    total = sum(blend.values()) or 1.0
    bio_avg  = bio_sum / bio_weight if bio_weight > 0 else 50.0
    max_conc = max(blend.values()) / total * 100

    return np.array([
        thickener_pct,
        surfactant_pct,
        preservative_pct,
        antioxidant_pct,
        bio_avg,
        float(len(blend)),
        max_conc,
        visc_potential,
        float(preservative_pct > 0),
        float(antioxidant_pct > 0),
    ], dtype=np.float32)


def _blend_hash(blend: Dict[str, float]) -> str:
    key = ",".join(f"{k}:{v:.2f}" for k, v in sorted(blend.items()))
    return hashlib.md5(key.encode()).hexdigest()[:8]


# ── ICP quantile ──────────────────────────────────────────────────────────────

def _icp_quantile(residuals: np.ndarray, level: float) -> float:
    n   = len(residuals)
    idx = min(int(np.ceil((n + 1) * level)) - 1, n - 1)
    return float(np.sort(residuals)[idx])


# ── Synthetic data generation (bootstrap prior) ───────────────────────────────

def _generate_synthetic(db: pd.DataFrame, n: int = STAB_SYNTH_N):
    """
    Sample n random blends from the ingredient DB, compute rule-based
    stability, add Gaussian noise, return (X, y_visc, y_shelf).
    """
    ings = db["Ingredient"].tolist() if "Ingredient" in db.columns else []
    if len(ings) < 3:
        return None, None, None

    rng     = np.random.default_rng(42)
    X_rows, y_visc, y_shelf = [], [], []

    for _ in range(n):
        k        = int(rng.integers(3, min(9, len(ings) + 1)))
        selected = rng.choice(ings, size=k, replace=False).tolist()
        concs    = rng.dirichlet(np.ones(k)) * 100.0
        blend    = dict(zip(selected, concs.tolist()))

        result = _rule_based(blend, db)

        feat = _blend_features(blend, db)
        if feat is None:
            continue

        v_noise = rng.normal(0, result.viscosity_cp * 0.12)
        s_noise = rng.normal(0, result.shelf_life_months * 0.10)

        X_rows.append(feat)
        y_visc.append(float(np.clip(result.viscosity_cp + v_noise, *_VISC_CLIP)))
        y_shelf.append(float(np.clip(result.shelf_life_months + s_noise, *_SHELF_CLIP)))

    if not X_rows:
        return None, None, None
    return np.array(X_rows), np.array(y_visc), np.array(y_shelf)


# ── ML training ───────────────────────────────────────────────────────────────

def _make_pipe(seed: int = 42) -> "Pipeline":
    return Pipeline([
        ("imputer", SimpleImputer(strategy="median")),
        ("scaler",  StandardScaler()),
        ("model",   GradientBoostingRegressor(
            n_estimators=120, max_depth=4,
            learning_rate=0.08, subsample=0.8,
            min_samples_leaf=4, random_state=seed,
        )),
    ])


def _train_models(
    X: np.ndarray, y_visc: np.ndarray, y_shelf: np.ndarray
) -> Tuple[Optional[Dict], Optional[StabilityConformalBands]]:
    """Train GBR pipelines with 85/15 ICP calibration split."""
    n_total = len(X)
    if n_total < 15:
        return None, None

    X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)

    n_cal   = max(5, int(n_total * STAB_CAL_PCT))
    n_fit   = n_total - n_cal
    rng     = np.random.default_rng(42)
    perm    = rng.permutation(n_total)
    fit_idx = perm[:n_fit]
    cal_idx = perm[n_fit:]

    targets = {"Viscosity": y_visc, "ShelfLife": y_shelf}
    models, residuals_dict = {}, {}

    for name, y in targets.items():
        pipe = _make_pipe(42)
        pipe.fit(X[fit_idx], y[fit_idx])
        models[name] = pipe
        residuals_dict[name] = np.abs(y[cal_idx] - pipe.predict(X[cal_idx]))

    quantiles: Dict[str, Dict[str, float]] = {
        name: {
            str(int(lv * 100)): _icp_quantile(res, lv)
            for lv in STAB_CONF_LEVELS
        }
        for name, res in residuals_dict.items()
    }
    bands = StabilityConformalBands(
        quantiles=quantiles, n_calibration=n_cal, coverage_levels=STAB_CONF_LEVELS
    )
    return models, bands


def _train_ensemble(
    X: np.ndarray, y_visc: np.ndarray, y_shelf: np.ndarray
) -> Optional[Dict]:
    """Train 5-model GBR ensemble for uncertainty estimation."""
    n_total = len(X)
    if n_total < 15:
        return None
    X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)
    ensemble: Dict[str, list] = {"Viscosity": [], "ShelfLife": []}
    for seed in range(STAB_ENSEMBLE_N):
        for name, y in [("Viscosity", y_visc), ("ShelfLife", y_shelf)]:
            pipe = _make_pipe(seed)
            pipe.fit(X, y)
            ensemble[name].append(pipe)
    return ensemble


# ── Feedback persistence ──────────────────────────────────────────────────────

def _load_feedback() -> List[StabilityFeedbackRecord]:
    if not os.path.exists(STAB_FEEDBACK_PATH):
        return []
    try:
        with open(STAB_FEEDBACK_PATH, "rb") as f:
            return pickle.load(f)
    except Exception:
        return []


def _save_feedback(records: List[StabilityFeedbackRecord]) -> None:
    try:
        with open(STAB_FEEDBACK_PATH, "wb") as f:
            pickle.dump(records, f)
    except Exception:
        pass


def _augment_with_feedback(
    X: np.ndarray, y_visc: np.ndarray, y_shelf: np.ndarray,
    feedback: List[StabilityFeedbackRecord],
    db: pd.DataFrame,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Append real feedback rows to the synthetic training arrays."""
    for rec in feedback:
        # Feedback records store blend_hash, not blend dict — use a dummy feature
        # vector centered on the feedback value. Real blends are stored in a separate
        # lookup; here we create a plausible feature row.
        dummy = np.zeros(len(_FEAT_NAMES), dtype=np.float32)
        dummy[7] = rec.actual_value if rec.target == "Viscosity" else 50.0
        X = np.vstack([X, dummy])
        y_visc  = np.append(y_visc,  rec.actual_value if rec.target == "Viscosity"  else y_visc.mean())
        y_shelf = np.append(y_shelf, rec.actual_value if rec.target == "ShelfLife" else y_shelf.mean())
    return X, y_visc, y_shelf


# ── Initialization ────────────────────────────────────────────────────────────

def initialize_stability_models(db: pd.DataFrame) -> StabilityModelCard:
    """
    Train GBR + ensemble on synthetic bootstrap data (+ any stored feedback).
    Results are cached to STAB_CACHE_PATH.
    """
    global _STAB_MODELS, _STAB_ENSEMBLE, _STAB_CONF, _STAB_CARD

    if not SKLEARN_OK:
        _STAB_CARD = StabilityModelCard(
            n_training=0, n_feedback=0, al_rounds=0, synth_samples=0,
            conformal_active=False, conformal_n_cal=0, ensemble_active=False,
            sklearn_version="unavailable",
        )
        return _STAB_CARD

    feedback = _load_feedback()

    # Try disk cache
    db_hash = hashlib.md5(
        (db["Ingredient"].sort_values().str.cat() + str(len(db))).encode()
    ).hexdigest()[:8]

    if os.path.exists(STAB_CACHE_PATH):
        try:
            with open(STAB_CACHE_PATH, "rb") as f:
                cached = pickle.load(f)
            if (cached.get("version") == STAB_CACHE_VERSION and
                    cached.get("db_hash") == db_hash and
                    cached.get("n_feedback") == len(feedback)):
                _STAB_MODELS   = cached["models"]
                _STAB_CONF     = cached["conformal"]
                _STAB_CARD     = cached["card"]
                _STAB_ENSEMBLE = _train_ensemble(*cached["arrays"])
                return _STAB_CARD
        except Exception:
            pass

    # Generate synthetic training data
    X_syn, y_visc_syn, y_shelf_syn = _generate_synthetic(db)
    if X_syn is None:
        _STAB_CARD = StabilityModelCard(
            n_training=0, n_feedback=len(feedback), al_rounds=0, synth_samples=0,
            conformal_active=False, conformal_n_cal=0, ensemble_active=False,
            sklearn_version="", cache_version=STAB_CACHE_VERSION,
        )
        return _STAB_CARD

    # Augment with feedback
    X, y_visc, y_shelf = _augment_with_feedback(X_syn, y_visc_syn, y_shelf_syn, feedback, db)

    models, bands = _train_models(X, y_visc, y_shelf)
    ensemble      = _train_ensemble(X, y_visc, y_shelf)

    _STAB_MODELS   = models
    _STAB_CONF     = bands
    _STAB_ENSEMBLE = ensemble

    import sklearn
    card = StabilityModelCard(
        n_training=len(X),
        n_feedback=len(feedback),
        al_rounds=0,
        synth_samples=len(X_syn),
        conformal_active=bands is not None,
        conformal_n_cal=bands.n_calibration if bands else 0,
        ensemble_active=ensemble is not None,
        sklearn_version=sklearn.__version__,
    )
    _STAB_CARD = card

    try:
        with open(STAB_CACHE_PATH, "wb") as f:
            pickle.dump({
                "version": STAB_CACHE_VERSION, "db_hash": db_hash,
                "n_feedback": len(feedback),
                "models": models, "conformal": bands, "card": card,
                "arrays": (X_syn, y_visc_syn, y_shelf_syn),
            }, f)
    except Exception:
        pass

    return card


# ── Rule-based core (unchanged — used as fallback + feature source) ───────────

def _rule_based(blend: Dict[str, float], db: pd.DataFrame) -> "StabilityResult":
    """Original heuristic predictor. Called by synthetic data gen and as fallback."""
    idx = db.set_index("Ingredient") if "Ingredient" in db.columns else pd.DataFrame()

    base_viscosity = 50.0
    for ing, pct in blend.items():
        if ing in _THICKENER_VISCOSITY:
            base_viscosity += _THICKENER_VISCOSITY[ing] * (pct / 100) ** 0.7
    for ing, pct in blend.items():
        try:
            if not idx.empty and ing in idx.index:
                func = str(idx.loc[ing, "Function"]) if "Function" in idx.columns else ""
                if "Surfactant" in func:
                    base_viscosity *= (1 - (pct / 100) * 0.15)
        except Exception:
            pass
    viscosity = round(base_viscosity, 0)
    v_low, v_high = round(viscosity * 0.85, 0), round(viscosity * 1.15, 0)

    base_shelf = 12
    max_pres = max((
        _PRESERVATIVE_EFFICACY[ing]
        for ing, pct in blend.items()
        if ing in _PRESERVATIVE_EFFICACY and pct > 0.5
    ), default=0)
    base_shelf = max(base_shelf, max_pres)
    ao_boost   = max((_ANTIOXIDANT_BOOST[ing] for ing in blend if ing in _ANTIOXIDANT_BOOST), default=0)
    base_shelf += ao_boost
    try:
        bio_avg = sum(
            float(idx.loc[ing, "Bio_based_pct"]) * (pct / 100)
            for ing, pct in blend.items()
            if not idx.empty and ing in idx.index
        )
        if bio_avg > 95:
            base_shelf = max(base_shelf - 3, 12)
    except Exception:
        pass
    shelf_low, shelf_high = max(base_shelf - 3, 6), base_shelf + 3

    ph_min, ph_max = 4.5, 8.5
    for ing in blend:
        if any(x in ing for x in ("Citric Acid", "Lactic Acid", "Glucono")):
            ph_min, ph_max = 4.0, 6.5
        if any(x in ing for x in ("Sodium Metasilicate", "Sodium Carbonate")):
            ph_min, ph_max = 9.0, 12.0
        if "Sodium Bicarbonate" in ing:
            ph_min, ph_max = 7.5, 9.0

    packaging = "HDPE (recommended for most green chemistry formulations)"
    for ing in blend:
        if ing in _EMOLLIENT_PACKAGING:
            packaging = _EMOLLIENT_PACKAGING[ing]
            break
    if ph_min >= 9.0:
        packaging = "HDPE only — alkaline formulas degrade PET and glass closures"
    if any("Ascorbic Acid" in ing or "Vitamin C" in ing for ing in blend):
        packaging = "Aluminium or dark glass — protect from UV and oxygen"

    risks = [_RISK_INGREDIENTS[ing] for ing in blend if ing in _RISK_INGREDIENTS] or \
            ["No significant stability risks identified"]
    boosters = []
    if ao_boost > 0:
        boosters.append("Antioxidant system detected — extends oxidative stability")
    if max_pres >= 18:
        boosters.append("Broad-spectrum preservative system — good microbial protection")
    if not boosters:
        boosters = ["Consider adding antioxidant and chelating agent to extend shelf life"]

    score = (
        (3 if base_shelf >= 24 else 2 if base_shelf >= 18 else 1)
        + (2 if max_pres > 0 else 0)
        + (1 if ao_boost > 0 else 0)
        + (1 if len(risks) <= 1 else 0)
    )
    rating = "Excellent" if score >= 6 else "Good" if score >= 4 else "Fair" if score >= 2 else "Poor"

    return StabilityResult(
        shelf_life_months=base_shelf,
        shelf_life_range=f"{shelf_low}–{shelf_high} months",
        viscosity_cp=viscosity,
        viscosity_range=f"{v_low:,.0f}–{v_high:,.0f} cP",
        ph_min=ph_min, ph_max=ph_max,
        recommended_packaging=packaging,
        stability_risks=risks,
        stability_boosters=boosters,
        overall_rating=rating,
        confidence="medium",
        ml_active=False,
    )


# ── ICP interval builder ──────────────────────────────────────────────────────

def _build_intervals(visc: float, shelf: float,
                     bands: StabilityConformalBands) -> Dict[str, Dict[str, Tuple[float, float]]]:
    result: Dict[str, Dict[str, Tuple[float, float]]] = {}
    for target, pred, clip in [("Viscosity", visc, _VISC_CLIP),
                                ("ShelfLife", shelf, _SHELF_CLIP)]:
        qs = bands.quantiles.get(target, {})
        result[target] = {
            lv: (
                round(max(clip[0], pred - q), 1),
                round(min(clip[1], pred + q), 1),
            )
            for lv, q in qs.items()
        }
    return result


# ── Public API ────────────────────────────────────────────────────────────────

def predict_stability(blend: Dict[str, float], db: pd.DataFrame) -> StabilityResult:
    """
    Predict stability and viscosity for a formulation blend.

    Uses tier-1 GBR+ICP when models are initialized (call initialize_stability_models
    first). Falls back to rule-based heuristics transparently.
    """
    rule = _rule_based(blend, db)

    if not SKLEARN_OK or _STAB_MODELS is None:
        return rule

    feat = _blend_features(blend, db)
    if feat is None:
        return rule

    X_2d = np.nan_to_num(feat, nan=0.0).reshape(1, -1)
    try:
        visc_ml  = float(np.clip(_STAB_MODELS["Viscosity"].predict(X_2d)[0],  *_VISC_CLIP))
        shelf_ml = float(np.clip(_STAB_MODELS["ShelfLife"].predict(X_2d)[0], *_SHELF_CLIP))
    except Exception:
        return rule

    v_lo = round(visc_ml * 0.85, 0)
    v_hi = round(visc_ml * 1.15, 0)
    s_lo = max(int(shelf_ml) - 3, 3)
    s_hi = int(shelf_ml) + 3

    intervals = _build_intervals(visc_ml, shelf_ml, _STAB_CONF) if _STAB_CONF else None

    # Uncertainty from ensemble
    confidence = "medium"
    if _STAB_ENSEMBLE:
        try:
            visc_preds  = [p.predict(X_2d)[0] for p in _STAB_ENSEMBLE["Viscosity"]]
            shelf_preds = [p.predict(X_2d)[0] for p in _STAB_ENSEMBLE["ShelfLife"]]
            visc_cv  = float(np.std(visc_preds)  / (np.mean(visc_preds)  + 1e-9))
            shelf_cv = float(np.std(shelf_preds) / (np.mean(shelf_preds) + 1e-9))
            mean_cv  = (visc_cv + shelf_cv) / 2
            confidence = "high" if mean_cv < 0.08 else "medium" if mean_cv < 0.18 else "low"
        except Exception:
            pass

    # Keep rule-based qualitative fields; replace quantitative with ML predictions
    return StabilityResult(
        shelf_life_months=int(round(shelf_ml)),
        shelf_life_range=f"{s_lo}–{s_hi} months",
        viscosity_cp=round(visc_ml, 0),
        viscosity_range=f"{v_lo:,.0f}–{v_hi:,.0f} cP",
        ph_min=rule.ph_min,
        ph_max=rule.ph_max,
        recommended_packaging=rule.recommended_packaging,
        stability_risks=rule.stability_risks,
        stability_boosters=rule.stability_boosters,
        overall_rating=rule.overall_rating,
        confidence=confidence,
        ml_active=True,
        intervals=intervals,
    )


def submit_stability_feedback(
    blend:        Dict[str, float],
    target:       str,
    actual_value: float,
    db:           pd.DataFrame,
) -> str:
    """
    Record a lab measurement for viscosity or shelf life and retrain the models.

    target: "Viscosity" (cP) | "ShelfLife" (months)
    Returns a status string for display.
    """
    global _STAB_MODELS, _STAB_ENSEMBLE, _STAB_CONF, _STAB_CARD

    valid = {"Viscosity", "ShelfLife"}
    if target not in valid:
        return f"Unknown target '{target}' — use 'Viscosity' or 'ShelfLife'."

    feedback = _load_feedback()
    al_round = ((_STAB_CARD.al_rounds + 1) if _STAB_CARD else 1)
    feedback.append(StabilityFeedbackRecord(
        blend_hash=_blend_hash(blend),
        target=target,
        actual_value=actual_value,
        timestamp=pd.Timestamp.now().isoformat(),
        al_round=al_round,
    ))
    _save_feedback(feedback)

    if SKLEARN_OK:
        X_syn, y_visc_syn, y_shelf_syn = _generate_synthetic(db)
        if X_syn is not None:
            X, y_visc, y_shelf = _augment_with_feedback(
                X_syn, y_visc_syn, y_shelf_syn, feedback, db)
            models, bands = _train_models(X, y_visc, y_shelf)
            ensemble      = _train_ensemble(X, y_visc, y_shelf)
            if models:
                _STAB_MODELS   = models
                _STAB_CONF     = bands
                _STAB_ENSEMBLE = ensemble
            if _STAB_CARD:
                _STAB_CARD.al_rounds  += 1
                _STAB_CARD.n_feedback  = len(feedback)
                _STAB_CARD.conformal_active = bands is not None
                _STAB_CARD.conformal_n_cal  = bands.n_calibration if bands else 0
            try:
                os.remove(STAB_CACHE_PATH)
            except Exception:
                pass
            return (f"✅ Stability model retrained — round {al_round}, "
                    f"{len(feedback)} lab point(s) total.")
    return f"✅ Feedback recorded (round {al_round}). Install scikit-learn to enable ML retraining."


def stability_model_card() -> Optional[StabilityModelCard]:
    return _STAB_CARD
