"""
计算分子空间电子密度
整个分子只有网格点，但是没有权重
对于每个原子，根据自身的坐标及权重插值到整个空间的坐标及权重
？ 波函数值是否也要插值呢？

空间中一个格点的权重是否应该对所有分子都一致？
"""
from pywfn.base import Mol
from pywfn.spaceprop import wfnfunc
from functools import cached_property,lru_cache
from pywfn.spaceprop import dftgrid
from pywfn.spaceprop import LineGrid,RectGrid,CubeGrid,MapsGrid
from pywfn import spaceprop
from pywfn import maths
import numpy as np

class Calculator(spaceprop.SpaceCaler):
    def __init__(self,mol:Mol) -> None:
        self.mol=mol
        self.wfnCaler=wfnfunc.Calculator(mol)
        self.PM=self.mol.PM.copy()
    
    def molDens_obt(self,grid:np.ndarray): #以分子轨道的方式计算电子密度
        """根据分子轨道电子密度计算分子电子密度"""
        obts=self.mol.O_obts
        dens=np.zeros(len(grid))
        wfns=self.wfnCaler.obtWfns(grid,obts) #分子轨道波函数
        for o,obt in enumerate(obts):
            wfn=wfns[o]
            dens+=wfn**2*self.mol.oE
        return dens
    
    def molDens_lib(self,grid:np.ndarray):
        """使用Fortran库计算电子密度"""
        from pywfn.maths import flib
        ngrid=len(grid)
        nmat=self.mol.CM.shape[0]
        nobt=len(self.mol.O_obts)
        obts=self.mol.O_obts
        CM=self.mol.CM[:,obts].copy()
        wfns=self.wfnCaler.atoWfns(grid) # 原子轨道波函数
        dens=flib.molDens(ngrid,nmat,nobt,CM,wfns)
        return dens*self.mol.oE
    
    def molDens_pro(self,grid:np.ndarray):
        """计算前体分子电子密度 promol"""
        pass
        
    def atmDens(self,grid:np.ndarray,atms:list[int]):
        """计算指定原子的电子密度加和，使用分子空间坐标"""
        nmat=self.mol.CM.shape[0]
        dens=np.zeros((len(grid)))
        for a,atm in enumerate(atms):
            atom=self.mol.atom(atm)
            u,l = atom.obtBorder
            atowfns=self.wfnCaler.atoWfns(grid)
            for i in range(u,l):
                wfn_i=atowfns[i]
                for j in range(nmat):
                    wfn_j=atowfns[j]
                    dens+=wfn_i*wfn_j*self.PM[i,j]
        return dens
    
    def vandSurf(self):
        p0,p1=self.mol.molBorder
        p0-=4
        p1+=4
        cube=CubeGrid().set(p0,p1,0.2)
        size,grid=cube.get()
        dens=self.molDens_lib(grid)
        verts,faces=self.isoSurf(dens,cube,0.1)
        return verts,faces