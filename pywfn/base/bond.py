from collections.abc import Iterator
from pywfn.base.atom import Atom
from pywfn import base
import numpy as np
# 如何获取一个键？
# Mol.getBond(a1,a2)
# Mol应该有bonds属性，是一个字典，索引是'a1-a2'和'a2-a1'都指向该键
from functools import cached_property
class Bond:
    def __init__(self,a1:Atom,a2:Atom) -> None:
        self.a1=a1
        self.a2=a2
        self._length:float=None
        self.idx:str=f'{a1.idx}-{a2.idx}'
        self.ats:list[int]=[a1.idx,a2.idx]
    
    @cached_property
    def length(self):
        """获取键长"""
        return np.linalg.norm(self.a2.coord-self.a1.coord)

    @cached_property
    def vector(self):
        """获取键向"""
        return self.a2.coord-self.a1.coord

    @cached_property
    def center(self):
        """获取键中心坐标"""
        return (self.a1.coord+self.a2.coord)/2
    
    def __repr__(self) -> str:
        return f'{self.a1.idx}-{self.a2.idx},{self.length:.4f}'
    
class Bonds:
    def __init__(self,mol:"base.Mol") -> None:
        self.mol=mol
        self.bonds:list[Bond]=[]
    
    def add(self,a1:Atom,a2:Atom):
        """添加一个键"""
        if a1.idx>=a2.idx:a1,a2=a2,a1
        bond=Bond(a1,a2)
        self.bonds.append(bond)
    
    def get(self,idx1:int,idx2:int):
        """获取一个键"""
        if idx1>=idx2:idx1,idx2=idx2,idx1
        for bond in self.bonds:
            if bond.a1.idx==idx1 and bond.a2.idx==idx2:
                return bond
        raise f'没有指定的键{idx1}-{idx2}'

    def __iter__(self) -> Iterator[Bond]:
        for bond in self.bonds:yield bond

    def __bool__(self)->bool:
        return len(self.bonds)!=0