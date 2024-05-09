"""
方向电荷差值的福井函数
"""
from pywfn.base import Mol
from pywfn.atomprop import dirProps
import numpy as np

class Calculator:
    def __init__(self,moln:Mol,mol0:Mol,molp:Mol):
        self.mols=[moln,mol0,molp]
        self.vects:list[np.ndarray]=None
        self.atoms:list[int]=None
        self.cals=[dirProps.Calculator(mol) for mol in self.mols]
        
    
    def calculate(self):
        for c,cal in enumerate(self.cals):
            chgs=cal.calculate()
            self.chgs[:,c]=chgs