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
        self.NM:np.ndarray=None
    
    def calculate(self)->np.ndarray:
        """
        返回与轨道系数矩阵对应的电子分布矩阵
        """
        CM=self.CM
        SM=self.mol.SM
        row,col=CM.shape
        oE=self.mol.oE
        NM=np.zeros(shape=(row,row,col))
        for a in range(row):
            for i in range(row):
                for j in self.mol.O_obts:
                    # 起始系数为0的可以直接跳过
                    NM[a,i,j]=CM[a,j]*CM[i,j]*SM[i,a]*oE
        # NM=NM.sum(axis=1)
        NM:np.ndarray=np.sum(NM,axis=1)[:,self.mol.O_obts] # 二维矩阵[n,occ]
        self.NM=NM
        NMo:np.ndarray=NM.sum(axis=0) # 每个轨道的电子数

        obtEngs=np.array(self.mol.obtEngs).reshape(1,-1)[:,self.mol.O_obts]
        # print(obtEngs)
        EM=NM*obtEngs
        return EM