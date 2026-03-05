"""
modules/chem_utils.py
RDKit helpers for IntelliForm v0.7.

Caching strategy:
  - get_mol_props() is called at DB load time, results stored in DataFrame
  - draw_mol() results are cached by SMILES string using st.cache_data
    so structures are only rendered once per session
"""
from functools import lru_cache
from typing import Dict, Optional

from rdkit import Chem
from rdkit.Chem import Draw, Descriptors
from PIL import Image


def get_mol_props(smiles: str) -> Dict[str, float]:
    """
    Compute molecular weight, LogP, and TPSA from a SMILES string.
    Returns zeros for invalid SMILES (never raises).
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return {
                'MW':   round(Descriptors.MolWt(mol),   1),
                'LogP': round(Descriptors.MolLogP(mol), 2),
                'TPSA': round(Descriptors.TPSA(mol),    1)
            }
    except Exception:
        pass
    return {'MW': 0.0, 'LogP': 0.0, 'TPSA': 0.0}


@lru_cache(maxsize=64)
def draw_mol(smiles: str, width: int = 200, height: int = 200) -> Optional[Image.Image]:
    """
    Render a molecule as a PIL Image.
    LRU-cached by SMILES+size so each structure is only drawn once.
    Returns None for invalid SMILES.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return Draw.MolToImage(mol, size=(width, height))
    except Exception:
        pass
    return None


def enrich_db(db):
    """
    Add MW, LogP, TPSA columns to the ingredients DataFrame in-place.
    Called once at startup.
    """
    for i, row in db.iterrows():
        props = get_mol_props(row['SMILES'])
        db.loc[i, 'MW']   = props['MW']
        db.loc[i, 'LogP'] = props['LogP']
        db.loc[i, 'TPSA'] = props['TPSA']
    return db
