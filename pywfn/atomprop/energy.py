"""
计算分子轨道内每个原子的能量
"""
from pywfn.base.mole import Mole
from pywfn.atomprop import lutils
from pywfn.utils import printer
from pywfn.maths.mol import projCM,engMat,piEleMat,eleMat
import numpy as np

class Calculator:
    def __init__(self,mole:Mole) -> None:
        self.mole=mole
        self.CM=self.mole.CM.copy()

    def atmEngs(self)->np.ndarray:
        """计算每个原子对应的轨道能量，原子电子能"""
        NM=eleMat(self.mole) # 电子分布矩阵
        EM=engMat(self.mole,NM) # 获取能量矩阵

        atoms=self.mole.atoms
        engs=np.zeros(len(atoms))
        for a,atom in enumerate(self.mole.atoms):
            u,l=atom.obtBorder
            atomEng=EM[u:l,:].sum()
            engs[a]=atomEng
        return engs
    
    def atmPiEngs(self):
        """计算投影后系数矩阵算出的能量，原子pi电子能"""
        NM=piEleMat(self.mole) # 获取能量矩阵
        EM=engMat(self.mole,NM)
        atoms=self.mole.atoms
        engs=np.zeros(len(atoms))
        for a,atom in enumerate(self.mole.atoms):
            u,l=atom.obtBorder
            atomEng=EM[u:l,:].sum()
            engs[a]=atomEng
        return engs