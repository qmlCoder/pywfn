"""
此脚本用来计算pi-Mulliken电子自旋
"""
import numpy as np
from pywfn.base import Mol,Atom
from pywfn.utils import printer
from pywfn.atomprop import lutils

class Calculator:
    def __init__(self,mol:"Mol") -> None:
        self.mol=mol
        self.vect:np.ndarray=None
        self.atom:int=None
    
    def get_Es(self,obts:list[int]): 
        """计算电子数量"""

        CM_=self.mol.projCM([self.atom],obts,[self.vect],zero=True)
        elects=lutils.get_ects(self.mol,obts,CM_)
        return elects

    def calculate(self):
        """计算所有原子的自旋"""
        if not self.mol.isOpenShell:
            printer.warn('非开壳层分子无法计算自旋')
            return
        # 首先要有alpha电子和beta电子对应的轨道系数
        obtNum=self.mol.CM.shape[0] # 系数矩阵行数，基函数数量
        a_obt=[i for i,e in enumerate(self.mol.obtEcts) if (e!=0 and i< obtNum)] # alpha 电子所在的轨道
        b_obt=[i for i,e in enumerate(self.mol.obtEcts) if (e!=0 and i>=obtNum)] # beta  电子所在的轨道
        aEs=np.array(self.get_Es(a_obt)) # alpha 电子数量
        bEs=np.array(self.get_Es(b_obt)) # beta  电子数量
        return (aEs-bEs)[self.atom-1]
        
    def printRes(self):
        resStr=self.resStr()
        printer.info('方向电子自旋分布: ')
        printer.res(resStr)
        
    def resStr(self):
        elect=self.calculate()
        return lutils.atomValueStr(self.mol,[self.atom],[elect])