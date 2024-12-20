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
from pywfn.data import radDens
import numpy as np

class Calculator(spaceprop.SpaceCaler):
    def __init__(self,mol:Mol) -> None:
        self.mol=mol
        self.wfnCaler=wfnfunc.Calculator(mol)
        self.PM=self.mol.PM.copy()
    
    def obtDens(self,grid:np.ndarray): # 以分子轨道的方式计算电子密度，可以计算每个分子轨道的电子密度
        """根据分子轨道电子密度计算分子电子密度"""
        obts=self.mol.O_obts
        dens=np.zeros((len(obts),len(grid)))
        wfns=self.wfnCaler.obtWfns(grid,obts) #分子轨道波函数
        for o,obt in enumerate(obts):
            wfn=wfns[o]
            dens[o]+=wfn**2*self.mol.oE
        return dens
    
    def molDens(self,grid:np.ndarray): # 这种方式计算分子的电子密度更快
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
    
    def proMolDens(self,grid:np.ndarray): # 第一种方式，使用插值计算
        """计算前体分子电子密度 promol"""
        dens=np.zeros((len(grid)))
        for atom in self.mol.atoms:
            radius=np.linalg.norm(grid-atom.coord,axis=1)
            dens+=radDens.get_radDens(atom.atomic,radius)
        return dens
    
    def proMolDens_v2(self,grid:np.ndarray,level=0): #第二种版本，使用公式计算
        """计算前体分子电子密度及梯度 promol"""
        dens0=np.zeros((len(grid)))
        dens1=np.zeros((len(grid),3))
        for atom in self.mol.atoms:
            rho0,rho1,dens2=radDens.get_radDens_v2(atom.atomic,grid,level)
            dens0+=rho0
            dens1+=rho1
        match level: # 根据level返回不同维度的密度
            case 0:
                return dens0
            case 1:
                return dens0,dens1
            case 2:
                return dens0,dens1,dens2
            case _:
                raise ValueError("level must be 0 or 1")
            
    def RDG(self,grid:np.ndarray,pro:bool=False): # reduced density gradient
        if pro: # 如果使用预分子
            dens0,dens1=self.proMolDens_v2(grid,1) # type: ignore
            return np.linalg.norm(dens1,axis=1)/(2*(3*np.pi**2)**(1/3)*dens0**(4/3))
        else: # 使用真实电子密度及梯度
            pass
    
    def IRI(self,grid:np.ndarray,pro:bool=False,a=1.1): # 
        if pro: # 如果使用预分子
            dens0,dens1=self.proMolDens_v2(grid,1) # type: ignore
            return np.linalg.norm(dens1,axis=1)/dens0**a
        else:
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
    
    def piMolDens(self,grid:np.ndarray):
        """计算分子空间π电子密度"""
        pass

    def vandSurf(self):
        p0,p1=self.mol.molBorder
        p0-=4
        p1+=4
        cube=CubeGrid().set(p0,p1,0.2)
        size,grid=cube.get()
        dens=self.molDens(grid)
        verts,faces=self.isoSurf(dens,cube,0.1)
        return verts,faces