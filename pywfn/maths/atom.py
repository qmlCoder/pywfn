import numpy as np

from pywfn.base.mole import Mole

counts={}

def get_sCont(mole:Mole,atm:int,obt:int):
    """获取某个原子轨道的贡献"""
    key=f'{mole.reader.path}-{atm}-{obt}'
    if key not in counts.keys():
        atom=mole.atom(atm)
        s=atom.OC[0,obt]**2
        As=np.sum(mole.CM**2,axis=0)[obt]
        counts[key]=s/As
    return counts[key]