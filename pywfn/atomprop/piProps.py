"""
该脚本计算指定原子的pi电子数
"""
from pywfn.base import Mol,Atom
from pywfn.maths import CM2PM
from pywfn.utils import printer
from pywfn.atomprop import lutils,dirProps,AtomCaler
from typing import Literal


class Calculator(AtomCaler):
    def __init__(self,mol:Mol) -> None:
        self.mol=mol
        self.caler=dirProps.Calculator(mol)
        self.chrg:Literal['mulliken','lowdin']=None
        self.prop:Literal['charge','spin']=None
        self.atoms:list[int]=None
        self.logTip:str=''

    def calculate(self):
        """
        计算所有原子的π电子
        如果指定方向，则计算方向电子
        每个原子都有应该根据自己的法向量方向来计算(当法向量都相同时相当于都一样)
        """
        atoms=self.mol.heavyAtoms
        atoms=lutils.atomIdxs(atoms)
        vects=lutils.get_vects(self.mol,atoms) # 每个原子都获得自己的法向量

        self.caler.chrg=self.chrg
        self.caler.prop=self.prop
        self.caler.vects=vects
        self.caler.atoms=atoms
        values=self.caler.calculate()
        return values

    def print(self,resStr:str):
        printer.info('所有非H原子的pi电子分布: ')
        printer.res(resStr)
    
    def resStr(self):
        elects=self.calculate()
        atoms=lutils.atomIdxs(self.mol.atoms)
        resStr=lutils.atomValueStr(self.mol,atoms,elects)
        resStr+=f'\ntotal: {sum(elects):.4f}'
        return resStr