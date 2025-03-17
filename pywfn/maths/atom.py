import numpy as np

from pywfn import base
from pywfn.base import Mole,Atom

def pola2cart(r,t,p):
    """极坐标转直角坐标，半径，仰角，转角"""
    z=r*np.sin(t)
    x=np.cos(t)*np.cos(p)
    y=np.cos(t)*np.sin(p)
    return x,y,z

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