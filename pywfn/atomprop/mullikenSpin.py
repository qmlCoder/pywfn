"""
此脚本用来计算mulliken自旋
"""
import numpy as np
from pywfn.base import Mol,Atom

from pywfn.utils import printer
from pywfn.atomprop import lutils

class Calculator:
    def __init__(self,mol:"Mol") -> None:
        self.mol=mol
    
    def get_Es(self,obts:list[int])->list[float]:
        elects=lutils.get_ects(self.mol,obts,self.mol.CM)
        return elects

    def calculate(self):
        """计算所有原子的自旋"""
        if not self.mol.isOpenShell:
            printer.warn('非开壳层分子无法计算自旋')
            return
        obtNum=self.mol.CM.shape[0] # 系数矩阵行数，基函数数量
        a_obt=[i for i,e in enumerate(self.mol.obtEcts) if (e!=0 and i< obtNum)]
        b_obt=[i for i,e in enumerate(self.mol.obtEcts) if (e!=0 and i>=obtNum)]
        a_Ects=np.array(self.get_Es(a_obt))
        b_Ects=np.array(self.get_Es(b_obt))
        return a_Ects-b_Ects
        
    def print(self,resStr:str):
        printer.info('mulliken 电子自旋分布:')
        printer.res(resStr)
        
    def resStr(self):
        elects=self.calculate()
        atoms=lutils.atomIdxs(self.mol.atoms)
        return lutils.atomValueStr(self.mol,atoms,elects)