"""计算一些性质的矩阵"""

from pywfn.base import Mol
from pywfn.maths import flib

import numpy as np

class Calculator:
    def __init__(self,mol:Mol) -> None:
        self.mol=mol
        self.CM=mol.CM.copy()
    
    def eleMat(self)->np.ndarray:
        """计算电子分布矩阵"""
        # 使用法向量可以计算每个分子的pi电子分布
        from pywfn.maths import flib
        
        SM=self.mol.SM
        oE=self.mol.oE
        
        obts=self.mol.O_obts
        nobt=len(obts)
        CM=self.mol.CM[:,obts].copy()
        nmat,nobt=CM.shape

        # NM=np.zeros(shape=(nmat,len(obts)))
        # for b in range(nmat):
        #     for j in range(nobt):
        #         v=0.0
        #         for i in range(nmat):
        #             # 起始系数为0的可以直接跳过
        #             v+=CM[b,j]*CM[i,j]*SM[i,b]*oE
        #         NM[b,j]=v
        # print(NM.sum())
        # return NM # 二维矩阵[n,occ]

        # print(CM[0,:])
        NM=flib.eleMat(nmat,nobt,CM,SM)*oE
        # print(NM.sum())
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