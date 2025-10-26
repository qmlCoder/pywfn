"""
定义存储分子结构信息的类 Geome
"""

import numpy as np
from pywfn import base
from pywfn.data.elements import elements
from pywfn import core

class Geome:
    def __init__(self) -> None:
        self.edited=False # 是否被编辑过
        self.atoms=base.Atoms(self)
        self.bonds=base.Bonds(self)
        self.mole:base.Mole|None=None # 绑定的分子对象
        self.core=core.base.Geome() # type: ignore

    def build(self,atms:list[int],xyzs:np.ndarray):
        """构建分子结构信息"""
        self.core.build(atms,xyzs)
    
    def addAtom(self,sym:str,xyz:np.ndarray):
        """添加一个原子"""
        self.atoms.add(sym,xyz)
        for atom in self.atoms:
            if atom.idx==len(self.atoms.atoms):continue
            r=np.linalg.norm(atom.coord-xyz)
            r1=elements[atom.symbol].radius
            r2=elements[sym].radius
            if r>(r1+r2)*1.1:continue
            self.bonds.add(atom.idx,len(self.atoms.atoms))
        self.edited=True

    @property
    def syms(self)->list[str]:
        return [atom.sym for atom in self.atoms]
    
    @property
    def xyzs(self):
        return np.array([atom.coord for atom in self.atoms])

    def __repr__(self) -> str:
        return f"{self.core}"
        
