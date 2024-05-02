"""
计算mayer键级
mayer键级计算公式为: 
$\sum_{a\in A}\sum_{b\in B}PS$
"""
import numpy as np
from pywfn.base import Atom,Mol
from pywfn.bondprop import Caler
from pywfn.maths import CM2PM
from pywfn.utils import printer

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
        PM=CM2PM(self.mol.CM,self.mol.O_obts,self.mol.oE)
        # 获取重叠矩阵
        SM=self.mol.SM
        PS=PM@SM
        OM=PS*PS.T
        a1,b1=atom1.obtBorder
        a2,b2=atom2.obtBorder
        order=np.sum(OM[a1:b1,a2:b2])
        return order
    
    def resStr(self) -> str:
        order=self.calculate()
        return f'{order:.4f}'
    
    def print(self):
        printer.res(self.resStr())

