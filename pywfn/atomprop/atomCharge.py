from pywfn.base import Mol,Atom
import numpy as np
from pywfn.data import Elements
from functools import lru_cache
elements=Elements()
from pywfn.utils import printer
from pywfn.atomprop import lutils,AtomCaler
from pywfn import maths
from typing import Literal

class Calculator(AtomCaler):
    def __init__(self,mol:"Mol"):
        self.logTip:str=''
        self.mol=mol
        self.chrg:Literal['mulliken','lowdin']='mulliken'
    
    def calculate(self)->np.ndarray:
        if self.chrg=='mulliken':
            return self.mulliken()
        if self.chrg=='lowdin':
            return self.lowdin()
    
    def mulliken(self):
        # 计算密度矩阵
        self.logTip='Mulliken电荷分布'
        PM=self.mol.PM
        # PS=PM*mol.SM
        # PSS=PS.sum(axis=0)
        # 矩阵乘法的迹的加和=矩阵对应元素乘积之和
        PS=PM@self.mol.SM
        EV=np.diagonal(PS)
        atoms=self.mol.atoms
        charges=np.zeros(len(atoms))
        for a,atom in enumerate(atoms):
            a1,a2=atom.obtBorder
            elect=EV[a1:a2].sum()
            charges[a]=atom.atomic-elect
        return charges

    def lowdin(self):
        """
        计算每个原子的lowdin电荷
        """
        self.logTip='Lowdin电荷分布'
        PM=self.mol.PM
        SM=self.mol.SM
        # 计算矩阵的1/2次方
        v, Q=np.linalg.eig(SM)
        V=np.diag(v)
        V_=V**0.5
        Q_=np.linalg.inv(Q)
        SM_half=Q@(V_@Q_)
        SPS=SM_half@(PM@SM_half)
        eleNums=np.diagonal(SPS)
        charges=np.zeros(len(self.mol.atoms))
        for a,atom in enumerate(self.mol.atoms):
            u,l=atom.obtBorder
            eleNum=eleNums[u:l].sum()
            charges[a]=atom.atomic-eleNum
        return charges
    
    def resStr(self)->str:
        """获取结果的打印内容"""
        satoms=lutils.atomIdxs(self.mol.atoms)
        charges=self.calculate()
        return lutils.atomValueStr(self.mol,satoms,charges)