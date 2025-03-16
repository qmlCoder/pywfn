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
from pywfn.spaceprop import LineGrid,RectGrid,CubeGrid,EarthGrid
from pywfn import spaceprop
from pywfn.data import radDens
from pywfn import utils
import numpy as np
import time

class Calculator(spaceprop.SpaceCaler):
    def __init__(self,mol:Mol) -> None:
        self.mol=mol
        self.wfnCaler=wfnfunc.Calculator(mol)
        self.PM=self.mol.PM.copy()
        
    
    def obtDens(self,grids:np.ndarray): # 以分子轨道的方式计算电子密度，可以计算每个分子轨道的电子密度
        """根据分子轨道电子密度计算分子电子密度"""
        obts=self.mol.O_obts
        dens=np.zeros((len(obts),len(grids)))
        wfns=self.wfnCaler.obtWfns(grids,obts) #分子轨道波函数
        for o,obt in enumerate(obts):
            wfn=wfns[o]
            dens[o]+=wfn**2*self.mol.oE
        return dens
    
    # def atmDens_(self,grids:np.ndarray,atms:list[int]|None=None):
    #     """计算指定原子的电子密度加和，使用分子空间坐标"""
    #     from pywfn.maths import flib
    #     if atms is None:atms=self.mol.atoms.atms
    #     nmat=self.mol.coefs.CM('car').shape[0]
    #     dens=np.zeros(shape=(len(atms),len(grids)))
    #     atowfns0,atowfns1,atowfns2=self.wfnCaler.atoWfns(grids,level=0) # 原子轨道波函数
    #     ngrid=grids.shape[0]
    #     atoDens=flib.atoDens(ngrid,nmat,self.mol.PM,atowfns0,atowfns1,atowfns2,0)[0]
    #     for a,atm in enumerate(atms):
    #         atom=self.mol.atom(atm)
    #         u,l = atom.obtBorder
    #         atoDens[u:l,:].sum(axis=0,out=dens[a,:])
    #     return dens
    
    def atmDens(self,grids:np.ndarray,atms:list[int]|None=None):
        """计算指定原子的电子密度加和，使用分子空间坐标"""
        from pywfn.maths import rlib
        if atms is None:atms=self.mol.atoms.atms
        dens=np.zeros(shape=(len(atms),len(grids)))
        xyzs,lmns,coes,alps=self.mol.basis.atoMap()
        atoDens=rlib.ato_rhos_rs(grids,xyzs,lmns,coes,alps,self.mol.PM,0)[0] # type: ignore
        atoDens=np.array(atoDens)
        for a,atm in enumerate(atms):
            atom=self.mol.atom(atm)
            u,l = atom.obtBorder
            dens[a,:]=atoDens[:,u:l].sum(axis=1)
        return dens
    
    # def molDens_(self,grids:np.ndarray,level:int)->tuple[np.ndarray,np.ndarray,np.ndarray]: # 这种方式计算分子的电子密度更快
    #     """使用Fortran库计算电子密度"""
    #     from pywfn.maths import flib
    #     ngrid=len(grids)
    #     obts=self.mol.O_obts
    #     CM=self.mol.coefs.CM('car')[:,obts].copy()
    #     nmat,nobt=CM.shape
    #     wfns0,wfns1,wfns2=self.wfnCaler.atoWfns(grids,level) # 原子轨道波函数
    #     dens0,dens1,dens2=flib.obtDens(ngrid,nmat,nobt,CM,wfns0,wfns1,wfns2,level)
    #     dens0=dens0.sum(axis=0)
    #     dens1=dens1.sum(axis=0)
    #     dens2=dens2.sum(axis=0)
    #     dens0*=self.mol.oE
    #     dens1*=self.mol.oE
    #     dens2*=self.mol.oE
    #     return dens0,dens1,dens2

    def molDens(self,grids:np.ndarray,level:int)->tuple[np.ndarray,np.ndarray,np.ndarray]:
        """使用Rust库计算电子密度

        Args:
            grids (np.ndarray): 要计算的格点
            level (int): 计算等级

        Returns:
            tuple[np.ndarray,np.ndarray,np.ndarray]: 电子密度、电子密度导数、电子密度Hessian
        """
        CM=self.mol.coefs.CM('car')[:,self.mol.O_obts].copy() # 占据轨道的系数矩阵
        xyzs,lmns,coes,alps=self.mol.basis.atoMap()
        dens0,dens1,dens2=rlib.mol_rhos_rs(grids.tolist(),xyzs,lmns,coes,alps,CM.tolist(),level) # type: ignore
        dens0=np.array(dens0)*self.mol.oE
        dens1=np.array(dens1)*self.mol.oE
        dens2=np.array(dens2)*self.mol.oE
        return dens0,dens1,dens2
    
    def proMolDens(self,grids:np.ndarray): # 第一种方式，使用插值计算
        """计算前体分子电子密度 promol"""
        dens=np.zeros((len(grids)))
        for atom in self.mol.atoms:
            radius=np.linalg.norm(grids-atom.coord,axis=1)
            dens+=radDens.get_radDens(atom.atomic,radius)
        return dens
    
    def proMolDens_v2(self,grid:np.ndarray,level=0): #第二种版本，使用公式计算
        """计算前体分子电子密度及梯度 promol"""
        dens0=np.zeros((len(grid),))
        dens1=np.zeros((len(grid),3))
        dens2=np.zeros((len(grid),3,3))
        for atom in self.mol.atoms:
            rho0,rho1,rho2=radDens.get_radDens_v2(atom.atomic,grid,level)
            dens0+=rho0
            dens1+=rho1
            dens2+=rho2
        return dens0,dens1,dens2
            
    def RDG(self,grids:np.ndarray,pro:bool=False): # reduced density gradient
        assert utils.chkArray(grids,[None,3]),"格点的形状应为: (n,3)"
        if pro: # 如果使用预分子
            dens0,dens1,_=self.proMolDens_v2(grids,1) # type: ignore
        else: # 使用真实电子密度及梯度
            dens0,dens1,_=self.molDens(grids,2)
        L=np.linalg.norm(dens1,axis=1)
        K=1/(2*(3*np.pi**2)**(1/3))
        rdg=K*L/dens0**(4/3)
        return rdg
    
    def IRI(self,grids:np.ndarray,pro:bool=False,a=1.1): # 
        assert utils.chkArray(grids,[None,3]),"格点的形状应为: (n,3)"
        if pro: # 如果使用预分子
            dens0,dens1,_=self.proMolDens_v2(grids,1) # type: ignore
        else:
            dens0,dens1,_=self.molDens(grids,1)
        L=np.linalg.norm(dens1,axis=1)
        # return L/dens0**a
        
        vals=L/dens0**a
        print(np.min(vals),np.max(vals))
        return vals
   
    
    
    def piMolDens(self,grid:np.ndarray):
        """计算分子空间π电子密度"""
        pass

    def vandSurf(self):
        """返回分子的范德华表面"""
        p0,p1=self.mol.molBorder
        shape,grids=spaceprop.CubeGrid().set_v1(p0,p1,0.3,4).get()
        vals=self.molDens(grids,0)[0]
        verts,faces=self.isoSurf(shape,grids,vals,0.04)
        assert verts is not None,"未找到范德华表面"
        assert faces is not None,"未找到范德华表面"
        return verts,faces