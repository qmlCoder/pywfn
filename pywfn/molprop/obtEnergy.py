"""
计算与分子轨道对应的电子分布
分子轨道中每个基函数的电子数量
"""
from pywfn.base import Mol
import numpy as np

class Calculator:
    def __init__(self,mol:Mol) -> None:
        self.mol=mol
        self.CM=mol.CM
    

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
                    NM[a,i,j]=CM[a,j]*CM[i,j]*SM[i,a]*oE
        NM=NM.sum(axis=1)
        sums=np.sum(NM,axis=0)
        NM=np.divide(NM,sums,out=np.zeros_like(NM),where=sums!=0) # 归一化
        obtEngs=np.array(self.mol.obtEngs)
        EM=NM*obtEngs
        return EM