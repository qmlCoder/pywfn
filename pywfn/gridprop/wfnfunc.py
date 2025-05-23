"""
计算分子轨道空间波函数
基函数的波函数
原子轨道的波函数
原子的波函数
"""
from pywfn.base import Mole
from pywfn import maths
import numpy as np
from pywfn.maths import cubeGrid
from pywfn.gridprop import lutils
from pywfn import gridprop

Array=np.ndarray

class Calculator(gridprop.SpaceCaler):
    def __init__(self,mol:Mole) -> None:
        self.mol=mol
        self.molPos=np.zeros((1,3)) # 初始坐标设为原点
        self.wfns:np.ndarray|None=None
        self.CM=self.mol.coefs.CM('car').copy()
        self.atms=self.mol.atoms.atms
    
    def obtWfns(self,grids:np.ndarray,obts:list[int])->np.ndarray:
        from pywfn.maths import rlib
        xyzs,lmns,coes,alps=self.mol.basis.atoMap()
        wfns0=[]
        for obt in obts:
            coefs=self.CM[:,obt].tolist() # 轨道系数
            # print('coes',coes)
            # print('alps',alps)
            wfn0,wfn1,wfn2=rlib.obt_wfns_rs(grids.tolist(),xyzs,lmns,coes,alps,coefs,0) # type: ignore
            # print(wfn0)
            wfns0.append(wfn0)
        return np.array(wfns0)

    
    # def obtWfns_(self,grid:np.ndarray,obts:list[int])->np.ndarray: #一次计算多个是最省性能的，而且多个也包含单个
    #     """
    #     计算分子轨道的波函数，为原子轨道的线性组合
    #     obt：分子轨道指标
    #     """
    #     atowfns,_,_=self.atoWfns(grid,level=0) # 所有原子轨道的波函数
    #     ngrid=grid.shape[0]
    #     nobt=len(obts)
    #     wfns=np.zeros(shape=(nobt,ngrid))
    #     for o,obt in enumerate(obts):
    #         coefs=self.CM[:,obt] # 轨道系数
    #         wfn=np.zeros(grid.shape[0])
    #         for c,coef in enumerate(coefs):
    #             wfn+=coef*atowfns[c]
    #         wfns[o]=wfn
    #     return wfns
    
    def atoWfns(self,grids:np.ndarray,level:int): # 所有原子轨道的波函数
        """计算所有原子轨道"""
        from pywfn.maths import rlib
        xyzs,lmns,coes,alps=self.mol.basis.atoMap()
        wfns0,wfns1,wfns2=rlib.ato_wfns_rs(grids,xyzs,lmns,coes,alps,level) # type: ignore
        wfns0=np.array(wfns0).swapaxes(0,1)
        wfns1=np.array(wfns1).swapaxes(0,1)
        wfns2=np.array(wfns2).swapaxes(0,1)
        return wfns0,wfns1,wfns2


    # def atoWfns_(self,grids:np.ndarray,level:int): # 所有原子轨道的波函数
    #     """计算所有原子轨道"""
    #     from pywfn.maths import flib
    #     ngrid=grids.shape[0]
    #     nmat=self.mol.coefs.CM('car').shape[0]
    #     # cords=self.mol.coords.copy()
    #     atms = self.mol.atoAtms
    #     shls = self.mol.atoShls
    #     syms = self.mol.atoSyms
    #     angs = self.mol.atoAngs
    #     basDict = self.mol.basis.dict
    #     coords=[]
    #     expl=[]
    #     coel=[]
    #     ncgs=[]
    #     cmax=0
    #     for i in range(nmat):
    #         atm=atms[i]
    #         atom=self.mol.atom(atm)
    #         # atomic=atom.atomic
    #         coords.append(atom.coord)
    #         shl=shls[i]
    #         ang=angs[i]
    #         key=f'{atom.idx}-{shl}-{ang}'
    #         assert key in basDict.keys(),"不存在的key"
    #         dat=np.array(basDict[key])
    #         exps=dat[:,0]
    #         coes=dat[:,1]
    #         expl.append(exps)
    #         coel.append(coes)
    #         if len(coes)>cmax:cmax=len(coes)
    #         ncgs.append(len(coes))
        
    #     coords=np.array(coords)
    #     expa=np.zeros(shape=(nmat,cmax))
    #     coea=np.zeros(shape=(nmat,cmax))
    #     for i,ncg in enumerate(ncgs):
    #         expa[i,:ncg]=expl[i]
    #         coea[i,:ncg]=coel[i]

    #     lmns = [self.mol.basis.sym2lmn(sym) for sym in syms]
    #     lmns=np.array(lmns)
    #     ncgs=np.array(ncgs)
        
    #     wfns0,wfns1,wfns2=flib.atoWfns(ngrid,grids,nmat,coords,cmax,ncgs,expa,coea,lmns,level)
    #     return wfns0,wfns1,wfns2

    def atmWfns(self,grid:np.ndarray,atms:list[int],obts:list[int])->np.ndarray: #一次计算多个是最省性能的，而且多个也包含单个
        """计算几个原子的波函数"""
        CM=self.mol.coefs.CM('car')
        nbas=CM.shape[0] # 基函数的数量
        obtAtms=self.mol.atoAtms
        atowfns,_,_=self.atoWfns(grid,level=0)
        nobt=len(obts)
        wfns=np.zeros(shape=(nobt,len(grid)))
        for o,obt in enumerate(obts):
            for i in range(nbas):
                if obtAtms[i] not in atms:continue
                coef=CM[i,obt]
                wfn=coef*atowfns[i]
                wfns[o]+=wfn
        return wfns