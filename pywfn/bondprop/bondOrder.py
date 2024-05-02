"""
键级也不止一种，都在这里实现吧
"""
from pywfn.base import Mol

import numpy as np

class Calculator:
    def __init__(self,mol:Mol) -> None:
        self.mol=mol

    def mayer(self):
        """计算mayer键级"""
        # 获取密度矩阵 P
        PM=self.mol.PM
        # 获取重叠矩阵
        SM=self.mol.SM
        PS=PM@SM
        OM=PS*PS.T
        orders=[]
        for bond in self.mol.bonds:
            u1,l1=bond.a1.obtBorder
            u2,l2=bond.a2.obtBorder
            order=np.sum(OM[u1:l1,u2:l2])
            orders.append([bond.atm1,bond.atm2,order])
        return orders
    
    def multiCenter(self,atms:list[int]):
        """计算多中心键级"""
        