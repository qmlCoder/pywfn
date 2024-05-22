"""
计算分子轨道空间波函数
基函数的波函数
原子轨道的波函数
原子的波函数
"""
from pywfn.base import Mol
from pywfn.maths import Gto
import numpy as np

from pywfn.data import sphGrid
weight = sphGrid.gridData[:, -1]
coords = sphGrid.gridData[:, :3]

class Calculator:
    def __init__(self,mol) -> None:
        self.mol:Mol=mol
        self.pos=np.zeros((1,3)) # 初始坐标设为原点
    
    def obtWfn(self,obt:int,pos:np.ndarray):
        """
        计算分子轨道的波函数，为原子轨道的线性组合
        obt：分子轨道指标
        coefs：线性组合系数
        """
        wfn=np.zeros(len(pos))
        coefs=self.mol.CM[:,obt]
        for c,coef in enumerate(coefs):
            wfn+=coef*self.atoWfn(c,pos)
        return wfn
    
    def atoWfn(self,i:int,pos:np.ndarray):
        """
        计算原子轨道的波函数，形成在轨道线性组合之前
        i:原子轨道指标
        """
        atms = self.mol.obtAtms
        shls = self.mol.obtShls
        syms = self.mol.obtSyms
        lmns = [self.mol.basis.sym2lmn(sym) for sym in syms]
        
        lmn = lmns[i]
        atm = atms[i]
        shl = shls[i]
        ang = sum(lmn)
        atmic = self.mol.atom(atm).atomic
        basis = self.mol.basis.get(atmic, shl, ang)
        exps = [b.exp for b in basis]
        coes = [b.coe for b in basis]
        pos_ = pos-self.mol.atom(atm).coord # 空间坐标-原子坐标=以原子为中心的空间坐标
        R2 = np.sum(pos_**2, axis=1)
        wfn = self.mol.gto.cgf(exps, coes, lmn, R2, pos_)  # 空间坐标-以原子为中心的坐标
        return wfn
