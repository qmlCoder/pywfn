"""
计算与分子轨道对应的电子分布
分子轨道中每个基函数的电子数量
"""
from pywfn.base import Mol
import numpy as np
from pywfn.utils import printer

class Calculator:
    def __init__(self,mol:Mol) -> None:
        self.mol=mol
        self.CM=mol.CM
        self.NM:np.ndarray|None=None
    
    def eleMat(self,CM:np.ndarray|None=None):
        """计算电子分布矩阵"""
        # 使用法向量可以计算每个分子的pi电子分布
        from pywfn.maths import flib
        if CM is None:CM=self.CM
        SM=self.mol.SM
        
        oE=self.mol.oE
        # row,col=CM.shape
        # NM=np.zeros(shape=(row,row,col))
        # for a in range(row):
        #     for i in range(row):
        #         for j in self.mol.O_obts:
        #             # 起始系数为0的可以直接跳过
        #             NM[a,i,j]=CM[a,j]*CM[i,j]*SM[i,a]*oE
        # NM:np.ndarray=np.sum(NM,axis=1)[:,self.mol.O_obts] # 二维矩阵[n,occ]
        nmat=CM.shape[0]
        obts=self.mol.O_obts
        nobt=len(obts)
        NM=flib.get_NM(nmat,nobt,CM[:,obts].copy(),SM)*oE
        return NM

    
    def engMat(self,CM:np.ndarray|None=None)->np.ndarray:
        """
        电子能量分布矩阵
        """
        if CM is None:CM=self.CM
        NM=self.eleMat(CM)
        # NMo:np.ndarray=NM.sum(axis=0) # 每个轨道的电子数
        # for o,obt in enumerate(self.mol.O_obts):
        #     print(obt+1,NMo[o])
        # print('总pi电子数量',np.sum(NMo))
        obtEngs=np.array(self.mol.obtEngs).reshape(1,-1)[:,self.mol.O_obts]
        # print(obtEngs)
        EM=NM*obtEngs
        return EM