"""
计算分子轨道内每个原子的能量
"""
from pywfn.base import Mol
from pywfn.molprop import obtEnergy
from pywfn.atomprop import lutils
from pywfn.utils import printer
import numpy as np

class Calculator:
    def __init__(self,mol:Mol) -> None:
        self.mol=mol
        self.caler=obtEnergy.Calculator(mol)
    
    def calculate(self)->np.ndarray:
        EM=self.caler.calculate() # 获取能量矩阵
        atoms=self.mol.atoms
        engs=np.zeros(len(atoms))
        for a,atom in enumerate(self.mol.atoms):
            u,l=atom.obtBorder
            atomEng=EM[u:l,:].sum()*self.mol.oE
            engs[a]=atomEng
        return engs
        
    def printRes(self):
        resStr=self.resStr()
        printer.info('原子能量分布: ')
        printer.res(resStr)
        
    def resStr(self):
        engs=self.calculate()
        atoms=[a.idx for a in self.mol.atoms]
        return lutils.atomValueStr(self.mol,atoms,engs)