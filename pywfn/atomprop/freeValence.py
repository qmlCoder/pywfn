"""
计算原子自由价
原子自由价=最大键级-该原子键级之和
这里的键级就先用DH
"""
from pywfn.bondprop import piDM
from pywfn.base import Mol,Atom
from pywfn.utils import printer
from pywfn.atomprop import AtomCaler
import numpy as np

STAND=1.6494416218465484

class Calculator(AtomCaler):
    def __init__(self,mol:Mol) -> None:
        self.mol=mol
        self.caler=piDM.Calculator(self.mol)
        self.vect:np.ndarray=None
        self.atom:int=None
        self.zero:bool=True
        self.keep:bool=True
        self.ins:bool=True
    
    def calculate(self):
        """计算一个原子的自由价，从1开始"""
        assert self.vect is not None,'未指定方向'
        atm=self.atom
        printer.info(f'计算原子 {atm} 自由价')
        centAtom=self.mol.atom(atm)
        valence1=STAND
        valence2=STAND
        self.caler.vect=self.vect.copy()
        self.caler.zero=self.zero
        self.caler.keep=self.keep
        self.caler.ins=self.ins

        for arouAtom in centAtom.neighbors:
            self.caler.bond=[atm,arouAtom.idx]
            orders=self.caler.calculate()
            order1,order2=orders
            bond=f'{centAtom.idx}-{arouAtom.idx}'
            printer.log(f'{bond} {order1:.6f} {order2:.6f}')
            valence1-=order1
            valence2-=order2
        if valence2==STAND:valence2=0.0
        return valence1,valence2
    
    
    def print(self,resStr:str):
        printer.info('原子自由价及其相关键级: ')
        printer.res(resStr)
    
    def resStr(self,idxs:list[int]):
        """结果的字符串形式"""
        valences=self.calculate(idxs)
        resStr=''
        for atom,(value1,value2) in zip(idxs,valences):
            resStr+=f'{atom:<8}{value1:<10.6f}{value2:<10.6f}\n'
        return resStr
