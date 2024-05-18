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
    
    def calculate(self,CM:np.ndarray=None)->np.ndarray:
        EM=self.caler.calculate(CM) # 获取能量矩阵
        atoms=self.mol.atoms
        engs=np.zeros(len(atoms))
        for a,atom in enumerate(self.mol.atoms):
            u,l=atom.obtBorder
            atomEng=EM[u:l,:].sum()
            engs[a]=atomEng
        return engs
    
    def dirEnergy(self):
        """计算投影后系数矩阵算出的能量"""
        from pywfn.atomprop import direction
        dirCaler=direction.Calculator(self.mol)
        atms=[]
        dirs=[]
        for atom in self.mol.atoms:
            normal=dirCaler.normal(atom.idx) # 原子的法向量
            if normal is None:continue
            atms.append(atom.idx)
            dirs.append(normal)
        CMp=self.mol.projCM(self.mol.O_obts,atms,dirs,True,False)
        engs=self.calculate(CMp)
        return engs

        
    def printRes(self):
        resStr=self.resStr()
        printer.info('原子能量分布: ')
        printer.res(resStr)
        
    def resStr(self):
        engs=self.calculate()
        atoms=[a.idx for a in self.mol.atoms]
        return lutils.atomValueStr(self.mol,atoms,engs)
    
    def onShell(self):
        engs=self.calculate()
        for i,eng in enumerate(engs):
            printer.res(f'{i}: {eng}')
        return