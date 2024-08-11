"""计算一些性质的矩阵"""

from pywfn.base import Mol
from pywfn.maths import flib

import numpy as np

class Calculator:
    def __init__(self,mol:Mol) -> None:
        self.mol=mol
        self.CM=mol.CM.copy()
    
    def eleMat(self,CM:np.ndarray|None=None)->np.ndarray:
        """计算电子分布矩阵"""
        # 使用法向量可以计算每个分子的pi电子分布
        from pywfn.maths import flib
        
        SM=self.mol.SM
        oE=self.mol.oE
        
        obts=self.mol.O_obts
        nobt=len(obts)
        if CM is None:
            CM=self.mol.CM[:,obts].copy()
        else:
            CM=CM[:,obts].copy()
        nmat,nobt=CM.shape

        NM=flib.eleMat(nmat,nobt,CM,SM)*oE
        # print(NM.sum())
        return NM
    
    def piEleMat(self):
        from pywfn.atomprop import direction
        dirCaler=direction.Calculator(self.mol)
        dirs=[]
        atms=[]
        for i in range(len(self.mol.atoms)):
            normal=dirCaler.normal(i+1)
            dirs.append(normal)
            atms.append(i+1)
        obts=self.mol.O_obts
        CMp=self.mol.projCM(obts,atms,dirs,False,False)
        NM=self.eleMat(CMp)
        return NM

    def engMat(self)->np.ndarray:
        """
        电子能量分布矩阵
        """
        if CM is None:CM=self.CM
        NM=self.eleMat(CM)
        obts=self.mol.O_obts
        obtEngs=np.array(self.mol.obtEngs).reshape(1,-1)[:,obts]
        EM=NM*obtEngs
        return EM