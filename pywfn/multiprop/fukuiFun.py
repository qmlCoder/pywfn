"""
根据三个分子的Mulliken电荷计算分子福井函数
"""
from pywfn.base import Mol
from pywfn.atomprop import charge
import numpy as np

class Calculator:
    def __init__(self,moln:Mol,mol0:Mol,molp:Mol) -> None:
        self.mols=[moln,mol0,molp]
        self.cals=[charge.Calculator(mol) for mol in self.mols]
        self.natm=len(mol0.atoms)
        self.chgs=np.zeros(shape=(self.natm,3)) # 记录所有原子的电荷
        self.dchg=np.zeros(shape=(self.natm,2)) # 记录电荷差值
    
    def calculate(self):
        for c,cal in enumerate(self.cals):
            chgs=cal.calculate()
            self.chgs[:,c]=chgs
        self.dchg[:,0]=-(self.chgs[:,0]-self.chgs[:,1])
        self.dchg[:,1]=-(self.chgs[:,1]-self.chgs[:,2])
        return self.dchg