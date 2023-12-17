import numpy as np

from pywfn.base import Mol

def molRange(mol:Mol)->tuple[np.ndarray,np.ndarray]:
    """分子在空间中的范围"""
    p0=mol.coords.min(axis=0)
    p1=mol.coords.max(axis=0)
    return p0,p1