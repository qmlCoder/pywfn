"""
计算分子空间电子密度
整个分子只有网格点，但是没有权重
对于每个原子，根据自身的坐标及权重插值到整个空间的坐标及权重
？ 波函数值是否也要插值呢？

空间中一个格点的权重是否应该对所有分子都一致？
"""
from pywfn.base import Mol
from pywfn.spaceProp import wfnfunc
from functools import cached_property,lru_cache
from pywfn.spaceProp import dftgrid,Line,Rect,Cube
from pywfn import maths
import numpy as np

class Calculator:
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
    
    def molDens_lib(self,grid):
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
        
    def atmDens(self,grid:np.ndarray,atms:list[int]):
        """计算指定原子的电子密度加和，使用分子空间坐标"""
        nmat=self.mol.CM.shape[0]
        dens=np.zeros((len(grid)))
        for a,atm in enumerate(atms):
            atom=self.mol.atom(atm)
            dens=np.zeros(len(grid))
            u,l = atom.obtBorder
            atowfns=self.wfnCaler.atoWfns(grid)
            for i in range(u,l):
                wfn_i=atowfns[i]
                for j in range(nmat):
                    wfn_j=atowfns[j]
                    dens+=wfn_i*wfn_j*self.PM[i,j]
        return dens
    
    def isoSurf(self): #获取电子密度等值面，顶点坐标以及面索引
        pass
    
    # def lineValue(self,line:Line,atms:list[int]):
    #     """获取一条直线上的数值"""
    #     size,grid=line.get()
    #     if not atms:
    #         dens=self.molDens_lib(grid)
    #     else:
    #         dens=self.atmDens(grid,atms)
    #     return dens.reshape(*size)
    
    # def rectValue(self,rect:Rect,atms:list[int]):
    #     """生成图片文件"""
    #     size,grid=rect.get()
    #     if not atms:
    #         dens=self.molDens_lib(grid)
    #     else:
    #         dens=self.atmDens(grid,atms)
    #     return dens.reshape(*size)
    
    # def cubeValue(self,cube:Cube,atms:list[int]):
    #     """生成图片文件"""
    #     size,grid=cube.get()
    #     if not atms:
    #         dens=self.molDens_lib(grid)
    #     else:
    #         dens=self.atmDens(grid,atms)
    #     return dens.reshape(*size)