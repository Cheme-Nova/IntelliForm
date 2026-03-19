"""
modules/chem_utils.py
RDKit helpers for IntelliForm.
"""
from functools import lru_cache
from typing import Dict, Optional

try:
    from rdkit import Chem
    from rdkit.Chem import Draw, Descriptors
    from PIL import Image
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False


def get_mol_props(smiles: str) -> Dict[str, float]:
    if not RDKIT_AVAILABLE:
        return {'MW': 0.0, 'LogP': 0.0, 'TPSA': 0.0}
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
def draw_mol(smiles: str, width: int = 200, height: int = 200):
    if not RDKIT_AVAILABLE:
        return None
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return Draw.MolToImage(mol, size=(width, height))
    except Exception:
        pass
    return None


def enrich_db(db):
    for i, row in db.iterrows():
        props = get_mol_props(row['SMILES'])
        db.loc[i, 'MW']   = props['MW']
        db.loc[i, 'LogP'] = props['LogP']
        db.loc[i, 'TPSA'] = props['TPSA']
    return db
