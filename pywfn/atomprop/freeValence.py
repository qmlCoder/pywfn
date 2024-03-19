"""
计算原子自由价
原子自由价=最大键级-该原子键级之和
这里的键级就先用DH
"""
from pywfn.bondorder import piDM
from pywfn.base import Mol,Atom
from pywfn.utils import printer
import numpy as np

STAND=1.6494416218465484

class Calculator:
    def __init__(self,mol:Mol) -> None:
        self.mol=mol
        self.caler=piDM.Calculator(self.mol)
        self.direct:np.ndarray=None
    
    def valence(self,idx:int):
        """计算一个原子的自由价，从1开始"""
        printer.info(f'计算原子 {idx} 自由价')
        assert self.direct is not None,'未指定方向'
        centAtom=self.mol.atom(idx)
        valence1=STAND
        valence2=STAND
        self.caler.direct=self.direct.copy()

        for arouAtom in centAtom.neighbors:
            self.caler.bond=[idx,arouAtom.idx]
            orders=self.caler.calculate()
            order1,order2=orders
            bond=f'{centAtom.idx}-{arouAtom.idx}'
            printer.log(f'{bond} {order1:.6f} {order2:.6f}')
            valence1-=order1
            valence2-=order2
        if valence2==STAND:valence2=0.0
        return valence1,valence2
    
    def calculate(self,idxs:list[int]):
        valences=[self.valence(idx) for idx in idxs]
        return valences
    
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
