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
        """计算每个原子对应的轨道能量"""
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
    
    def onShell(self):
        engs=self.calculate()
        while True:
            printer.options('原子能量',{
                '1':'原子轨道能',
                '2':'方向原子能'
            })
            opt=input('输入要计算的原子能:')
            if opt=='1':
                energy=self.calculate()
                for i,eng in enumerate(energy):
                    printer.res(f'{i+1}: {eng}')
            elif opt=='2':
                energy=self.dirEnergy()
                for i,eng in enumerate(engs):
                    printer.res(f'{i}: {eng}')
            else:
                break