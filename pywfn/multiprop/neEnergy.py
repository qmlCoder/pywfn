"""
计算亲核亲电能
"""
from pywfn.base import Mol
from pywfn.atomprop import atomEnergy

import numpy as np

class Calculator:
    def __init__(self,moln:Mol,mol0:Mol,molp:Mol) -> None:
        self.mols=[moln,mol0,molp]
        self.cals=[atomEnergy.Calculator(m) for m in self.mols]
        self.natm=len(mol0.atoms)
        self.engs=np.zeros((self.natm,3))
        self.deng=np.zeros((self.natm,2))
    
    def calculate(self):
        for c,cal in enumerate(self.cals):
            engs=cal.calculate()
            self.engs[:,c]=engs
        self.deng[:,0]=self.engs[:,0]-self.engs[:,1]
        self.deng[:,1]=self.engs[:,1]-self.engs[:,2]
        return self.deng