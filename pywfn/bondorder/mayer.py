"""
计算mayer键级
mayer键级计算公式为: 
$\sum_{a\in A}\sum_{b\in B}PS$
"""
import numpy as np
from pywfn.base import Atom,Mol
from pywfn.bondorder import Caler

class Calculator(Caler):
    def __init__(self,mol:"Mol"):
        self.mol=mol
        self.bond:list[int]=None

    def calculate(self)->float:
        idx1,idx2=self.bond
        atom1=self.mol.atom(idx1)
        atom2=self.mol.atom(idx2)
        """计算两原子之间的mayer键级"""
        # 获取密度矩阵 P
        PM=self.mol.PM
        # 获取重叠矩阵
        SM=self.mol.SM
        PS=PM@SM

        a_1,a_2=atom1.obtRange
        b_1,b_2=atom2.obtRange

        order=np.sum(PS[a_1:a_2,b_1:b_2]*PS[b_1:b_2,a_1:a_2].T)
        return order
