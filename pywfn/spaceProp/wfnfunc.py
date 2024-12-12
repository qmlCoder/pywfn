"""
计算分子轨道空间波函数
基函数的波函数
原子轨道的波函数
原子的波函数
"""
from pywfn.base import Mol
from pywfn import maths
import numpy as np
from pywfn.maths import flib
from pywfn.maths import cubeGrid
from pywfn.spaceprop import lutils
from pywfn import spaceprop
from pywfn.maths import march
Array=np.ndarray

class Calculator(spaceprop.SpaceCaler):
    def __init__(self,mol:Mol) -> None:
        self.mol=mol
        self.molPos=np.zeros((1,3)) # 初始坐标设为原点
        self.wfns:np.ndarray|None=None
        self.CM=self.mol.CM.copy()
        self.atms=self.mol.atoms.atms
    
    def obtWfns(self,grid:np.ndarray,obts:list[int])->np.ndarray: #一次计算多个是最省性能的，而且多个也包含单个
        """
        计算分子轨道的波函数，为原子轨道的线性组合
        obt：分子轨道指标
        """
        atowfns=self.atoWfns(grid) #所有原子轨道的波函数
        ngrid=grid.shape[0]
        nobt=len(obts)
        wfns=np.zeros(shape=(nobt,ngrid))
        for o,obt in enumerate(obts):
            coefs=self.CM[:,obt] # 轨道系数
            wfn=np.zeros(grid.shape[0])
            for c,coef in enumerate(coefs):
                wfn+=coef*atowfns[c]
            wfns[o]=wfn
        return wfns
    
    def atoWfns(self,grid:np.ndarray): #所有原子轨道的波函数
        """计算所有原子轨道"""
        ngrid=grid.shape[0]
        nmat=self.mol.CM.shape[0]
        # cords=self.mol.coords.copy()
        atms = self.mol.obtAtms
        shls = self.mol.obtShls
        syms = self.mol.obtSyms
        angs = self.mol.obtAngs
        basDict = self.mol.basis.dict
        coords=[]
        expl=[]
        coel=[]
        ncgs=[]
        cmax=0
        for i in range(nmat):
            atm=atms[i]
            atom=self.mol.atom(atm)
            atomic=atom.atomic
            coords.append(atom.coord)
            shl=shls[i]
            ang=angs[i]
            key=f'{atomic}-{shl}-{ang}'
            dat=np.array(basDict[key])
            exps=dat[:,0]
            coes=dat[:,1]
            expl.append(exps)
            coel.append(coes)
            if len(coes)>cmax:cmax=len(coes)
            ncgs.append(len(coes))
        
        coords=np.array(coords)
        expa=np.zeros(shape=(nmat,cmax))
        coea=np.zeros(shape=(nmat,cmax))
        for i,ncg in enumerate(ncgs):
            expa[i,:ncg]=expl[i]
            coea[i,:ncg]=coel[i]

        lmns = [self.mol.basis.sym2lmn(sym) for sym in syms]
        lmns=np.array(lmns)
        ncgs=np.array(ncgs)
        
        wfns=flib.atoWfns(ngrid,grid,nmat,coords,cmax,ncgs,expa,coea,lmns)
        return wfns

    def atmWfns(self,grid:np.ndarray,atms:list[int],obts:list[int])->np.ndarray: #一次计算多个是最省性能的，而且多个也包含单个
        """计算几个原子的波函数"""
        
        nbas=self.mol.CM.shape[0] # 基函数的数量
        obtAtms=self.mol.obtAtms
        atowfns=self.atoWfns(grid)
        nobt=len(obts)
        wfns=np.zeros(shape=(nobt,len(grid)))
        for o,obt in enumerate(obts):
            for i in range(nbas):
                if obtAtms[i] not in atms:continue
                coef=self.mol.CM[i,obt]
                wfn=coef*atowfns[i]
                wfns[o]+=wfn
        return wfns
    
def meshgrid(xr,yr):
    xs,ys=np.meshgrid(xr,yr)
    xs=xs.flatten()
    ys=ys.flatten()
    return xs,ys