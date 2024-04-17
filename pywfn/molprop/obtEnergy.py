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
                    NM[a,i,j]=CM[a,j]*CM[i,j]*SM[i,a]*oE
        # NM=NM.sum(axis=1)
        NM:np.ndarray=np.sum(NM,axis=1)[:,self.mol.O_obts] # 二维矩阵[n,occ]
        self.NM=NM
        NMo:np.ndarray=NM.sum(axis=0) # 每个轨道的电子数
        # print('-'*50)
        # print('分子总电子数',NMo.sum()) # 分子的总电子数
        # assert np.min(NM)>=-1e-3,'电子数不能小于0'
        # NM=np.divide(NM,NMo,out=np.zeros_like(NM),where=NMo!=0) # 归一化
        obtEngs=np.array(self.mol.obtEngs).reshape(1,-1)[:,self.mol.O_obts]
        # print(obtEngs)
        EM=NM*obtEngs
        # print('总电子数',np.sum(NM))
        # for a,atom in enumerate(self.mol.atoms):
        #     u,l=atom.obtBorder
        #     print(f'原子{a+1}电子数',NM[u:l,:].sum())
            # print(f'原子{a+1}电子能',EM[u:l,:].sum())
        return EM