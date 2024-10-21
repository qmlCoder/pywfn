"""
计算分子轨道内每个原子的能量
"""
from pywfn.base import Mol
from pywfn.atomprop import lutils
from pywfn.utils import printer
from pywfn.maths.mol import projCM,engMat,piEleMat,eleMat
import numpy as np

class Calculator:
    def __init__(self,mol:Mol) -> None:
        self.mol=mol
        self.CM=self.mol.CM.copy()
    
    def atmEngs(self)->np.ndarray:
        """计算每个原子对应的轨道能量，原子电子能"""
        NM=eleMat(self.mol) # 电子分布矩阵
        EM=engMat(self.mol,NM) # 获取能量矩阵
        
        atoms=self.mol.atoms
        engs=np.zeros(len(atoms))
        for a,atom in enumerate(self.mol.atoms):
            u,l=atom.obtBorder
            atomEng=EM[u:l,:].sum()
            engs[a]=atomEng
        return engs
    
    def atmPiEngs(self):
        """计算投影后系数矩阵算出的能量，原子pi电子能"""
        NM=piEleMat(self.mol) # 获取能量矩阵
        EM=engMat(self.mol,NM)
        atoms=self.mol.atoms
        engs=np.zeros(len(atoms))
        for a,atom in enumerate(self.mol.atoms):
            u,l=atom.obtBorder
            atomEng=EM[u:l,:].sum()
            engs[a]=atomEng
        return engs
    
    def onShell(self):
        while True:
            printer.options('原子能量',{
                '1':'原子轨道能',
                '2':'方向原子能'
            })
            opt=input('输入要计算的原子能:')
            if opt=='1':
                energy=self.atmEngs()
                for i,eng in enumerate(energy):
                    printer.res(f'{i+1}: {eng}')
            elif opt=='2':
                engs=self.atmPiEngs()
                for i,eng in enumerate(engs):
                    printer.res(f'{i}: {eng}')
            else:
                break