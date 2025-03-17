import numpy as np

from pywfn.base import Mole

counts={}

def get_sCont(mol:Mole,atm:int,obt:int):
    """获取某个原子轨道的贡献"""
    key=f'{mol.reader.path}-{atm}-{obt}'
    if key not in counts.keys():
        atom=mol.atom(atm)
        s=atom.OC[0,obt]**2
        As=np.sum(mol.CM**2,axis=0)[obt]
        counts[key]=s/As
    return counts[key]