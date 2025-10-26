

from collections.abc import Iterator
from pywfn import base
import numpy as np
from functools import cached_property

class Bond:
    """
    键对象
    """
    def __init__(self,geome:"base.Geome",atm1:int,atm2:int) -> None:
        self.geome=geome
        self.atm1=atm1
        self.atm2=atm2
        self._length:float|None=None
        self.key:str=f'{atm1}-{atm2}'
        self.ats:tuple[int,int]=(atm1,atm2)
    
    @property
    def mole(self):
        assert self.geome.mole is not None,"未绑定分子"
        return self.geome.mole
    
    @property
    def a1(self):
        return self.mole.atom(self.atm1)

    @property
    def a2(self):
        return self.mole.atom(self.atm2)

    @cached_property
    def length(self):
        """获取键长"""
        length=np.linalg.norm(self.a2.coord-self.a1.coord)
        return float(length)

    @cached_property
    def vector(self):
        """获取键向"""
        return self.a2.coord-self.a1.coord

    @cached_property
    def center(self):
        """获取键中心坐标"""
        return (self.a1.coord+self.a2.coord)/2
    
    @property
    def btype(self)->int:
        """获取键类型
        1:单键
        2:双键
        3:三键
        4:芳香
        """
        blen=self.length
        if blen<1.38*1.889:
            return 1
        else:
            return 2
    
    def __repr__(self) -> str:
        return f'{self.a1.idx:>3}-{self.a2.idx:>3},{self.length:>10.6f}'
    
class Bonds:
    def __init__(self,geome:"base.Geome") -> None:
        self.geome=geome
        self.bonds:list[Bond]=[]
    
    def add(self,idx1:int,idx2:int):
        """添加一个键"""
        if idx1>=idx2:idx1,idx2=idx2,idx1
        bond=Bond(self.geome,idx1,idx2)
        self.bonds.append(bond)
        return bond
    
    def get(self,idx1:int,idx2:int):
        """获取一个键"""
        if idx1>=idx2:idx1,idx2=idx2,idx1
        for bond in self.bonds:
            if bond.a1.idx==idx1 and bond.a2.idx==idx2:
                return bond
        return None
    
    def pop(self,idx1:int,idx2:int):
        """删除一个键"""
        for bond in self.bonds:
            if bond.ats==[idx1,idx2]:
                self.bonds.remove(bond)
                return bond
        return None # 没有删除成功
    
    
    @property
    def keys(self):
        return [bond.key for bond in self.bonds]
    
    @property
    def ats(self):
        return [bond.ats for bond in self.bonds]
    
    @property
    def num(self):
        return len(self.bonds)

    def __iter__(self) -> Iterator[Bond]:
        for bond in self.bonds:yield bond

    def __bool__(self)->bool:
        return len(self.bonds)!=0