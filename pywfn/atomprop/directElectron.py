from pywfn.base import Mol,Atom
from pywfn.atomprop import lutils
from pywfn.utils import printer
import numpy as np

class Calculator:
    def __init__(self,mol:Mol) -> None:
        self.mol=mol
        self.vect:np.ndarray=None # 指定的方向
        self.atom:int=None # 选择的原子(从1开始)


    def calculate(self):
        obts=self.mol.O_obts
        atoms=[self.atom]
        vects=[self.vect]
        CM_=self.mol.projCM(atoms,obts,vects,zero=True)
        ects=lutils.get_ects(self.mol,obts,CM_)
        return ects[self.atom-1]

    def printRes(self):
        resStr=self.resStr()
        printer.info(f'原子{self.atom}的方向电子数: ')
        printer.res(resStr)
    
    def resStr(self)->str:
        return f'{self.calculate():.4f}'