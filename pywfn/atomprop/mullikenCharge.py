from pywfn.base import Mol,Atom
import numpy as np
from pywfn.data import Elements
from functools import lru_cache
elements=Elements()
from pywfn.utils import printer
from pywfn.atomprop import lutils
class Calculator:
    def __init__(self,mol:"Mol"):
        self.mol=mol
    
    @lru_cache
    def calculate(self)->list[float]:
        elects=lutils.get_ects(self.mol,self.mol.O_obts,self.mol.CM)
        charges=[atom.atomic-elect for atom,elect in zip(self.mol.atoms,elects)]
        return charges
            
    def print(self,result):
        printer.info('Mulliken 电荷分布')
        printer.res(result)
    
    def resStr(self)->str:
        """获取结果的打印内容"""
        satoms=lutils.atomIdxs(self.mol.atoms)
        charges=self.calculate()
        return lutils.atomValueStr(self.mol,satoms,charges)