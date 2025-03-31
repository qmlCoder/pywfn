"""
定义存储分子结构信息的类 Geome
"""

import numpy as np

from pywfn import base
from pywfn.base.atom import Atoms
from pywfn.base.bond import Bonds
from pywfn.utils import chkArray
from pywfn.data.elements import elements

class Geome:
    def __init__(self,syms:list[str],xyzs:np.ndarray) -> None:
        self.syms:list[str]=syms # 原子符号
        self.xyzs:np.ndarray=xyzs # 原子坐标
        chkArray(xyzs,[None,3])
        self.edited=False # 是否被编辑过
        self.atoms=Atoms(self)
        self.bonds=Bonds(self)
        self.mol:"base.Mole|None"=None # 绑定的分子对象
    
    def update(self,syms,xyzs):
        """更新分子结构信息"""
        self.syms=syms
        self.xyzs=xyzs
        chkArray(xyzs,[None,3])
        self.edited=True # 被编辑过
        self.build()
    
    def build(self):
        """构建分子结构信息"""
        self.atoms.atoms.clear()
        self.bonds.bonds.clear()
        for sym,xyz in zip(self.syms,self.xyzs):
            self.atoms.add(sym,xyz)
        for atom1 in self.atoms:
            for atom2 in self.atoms:
                if atom1.idx>=atom2.idx:continue
                r=np.linalg.norm(atom2.coord-atom1.coord)
                r1=elements[atom1.symbol].radius
                r2=elements[atom2.symbol].radius
                if r>(r1+r2)*1.1:continue
                self.bonds.add(atom1.idx,atom2.idx)
        return self
    
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

    def __repr__(self) -> str:
        text=''
        for atom in self.atoms:
            sym=atom.sym
            x,y,z=atom.coord
            text+=f'{sym:>3}{x:>10.4f}{y:>10.4f}{z:>10.4f}\n'
        return text[:-1]
        
