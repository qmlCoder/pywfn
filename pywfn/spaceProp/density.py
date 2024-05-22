"""
计算分子空间电子密度
整个分子只有网格点，但是没有权重
对于每个原子，根据自身的坐标及权重插值到整个空间的坐标及权重
？ 波函数值是否也要插值呢？
"""
from pywfn.base import Mol
from pywfn.spaceProp import wfnfunc
from pywfn.data import sphGrid
from scipy.interpolate import griddata

import numpy as np

class Calculator:
    def __init__(self,mol:Mol) -> None:
        self.mol=mol
        self.wfnCaler=wfnfunc.Calculator(mol)
        self.points=sphGrid.gridData[:,:-1] # 网格的坐标
        self.weights=sphGrid.gridData[:,-1] # 网格的权重
    
    def molPos(self):
        """整个分子的网格点坐标"""
        natm=len(self.mol.atoms)
        molPos=[self.atmPos(i+1) for i in range(natm)]
        return np.vstack(molPos)

    def atmPos(self,atm:int):
        """单个原子的网格点坐标"""
        return self.points.copy()+self.mol.atom(atm).coord
    
    def atmWeight(self,atm:int):
        """单个原子在分子网格中的权重"""
        atmPos=self.atmPos(atm)
        molPos=self.molPos()
        # 根据插值获取权重
        atmWeight=griddata(atmPos,self.weights,molPos,method='linear')
        return atmWeight

    def molDens_obt(self,pos:np.ndarray):
        """计算根据分子轨道电子密度计算的分子电子密度"""
        obts = self.mol.O_obts
        dens=np.zeros(len(pos))
        for o in obts:
            wfn=self.wfnCaler.obtWfn(o,pos)
            dens+=wfn**2*self.mol.oE
        return dens
    
    def molDens_atm(self,pos:np.ndarray):
        """根据原子电子密度计算的分子密度"""
        dens=np.zeros(len(pos))
        for atom in self.mol.atoms:
            dens+=self.atmDens(atom.idx,pos)
        return dens
    
    def atmDens(self,atm:int,pos:np.ndarray):
        """计算原子电子密度"""
        nmat=self.mol.CM.shape[0]
        atom=self.mol.atom(atm)
        dens=np.zeros(len(pos))
        u,l = atom.obtBorder
        for i in range(nmat):
            wfn_i=self.wfnCaler.atoWfn(i,pos)
            for j in range(u,l):
                wfn_j=self.wfnCaler.atoWfn(j,pos)
                dens+=wfn_i*wfn_j*self.mol.PM[i,j]
        return dens