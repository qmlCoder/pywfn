"""
计算分子的芳香性
使用pi键级的标准差表示
"""
from pywfn.base import Mole
from pywfn.atomprop import charge
from pywfn.bondprop import order as orderProp
import numpy as np
from pywfn.cli import Shell
from pywfn.utils import printer

class Calculator:
    def __init__(self, mol:Mole):
        self.mol = mol
        self.ratio=0.5
    
    def pisd(self,ring:list[int]|None=None): # 版本1 直接用键级标准差
        caler=orderProp.Calculator(self.mol)
        bonds,orders=caler.pi_pocv()
        # orders=result[:,-1]
        forders=[] # 过滤掉C-H键
        for i,((a1,a2),order) in enumerate(zip(bonds,orders)):
            if self.mol.atom(a1).atomic==1:continue
            if self.mol.atom(a2).atomic==1:continue
            if len(self.mol.atom(a1).neighbors)==4:continue
            if len(self.mol.atom(a2).neighbors)==4:continue
            if ring is None:
                print(f'{a1:>2}-{a2:>2}:{order:>10.4f}')
                forders.append(order)
            elif (a1 in ring and a2 in ring):
                print(f'{a1:>2}-{a2:>2}:{order:>10.4f}')
                forders.append(order)
            
        if len(forders)==0:
            printer.warn('没有键')
            return 0
        return np.std(forders).item()
    
    def pimsd(self,ring:list[int]|None=None,ratio=0.5): # 版本2 使用键级均值和标准差
        caler=orderProp.Calculator(self.mol)
        bonds,orders=caler.pi_pocv()
        forders=[] # 过滤掉C-H键
        for i,((a1,a2),order) in enumerate(zip(bonds,orders)):
            if self.mol.atom(a1).atomic==1:continue
            if self.mol.atom(a2).atomic==1:continue
            if ring is None:
                print(f'{a1:>2}-{a2:>2}:{order:>10.4f}')
                forders.append(order)
            elif (a1 in ring and a2 in ring):
                print(f'{a1:>2}-{a2:>2}:{order:>10.4f}')
                forders.append(order)
        mean=np.mean(forders)
        stds=np.std(forders)
        return (ratio*mean-(1-ratio)*stds).item()

    def pimed(self): # 使用键级类比于HOMED方法
        D=0.2
        caler=orderProp.Calculator(self.mol)
        bonds,orders=caler.pi_pocv()
        forders=[] # 过滤掉C-H键
        for i,((a1,a2),order) in enumerate(zip(bonds,orders)):
            if self.mol.atom(a1).atomic==1:continue
            if self.mol.atom(a2).atomic==1:continue
            forders.append((order-0.66)**2)
        nbond=self.mol.bonds.num
        result=1-sum(forders)/(nbond*D)
        return result
    
    def homed(self):
        D=0.2
        idea=1.39645
        vals=[]
        for bond in self.mol.bonds:
            atm1,atm2=bond.atm1,bond.atm2
            if self.mol.atom(atm1).atomic==1:continue
            if self.mol.atom(atm2).atomic==1:continue
            vals.append((bond.length/1.889-idea)**2)
        val=sum(vals)
        return 1-val/(self.mol.bonds.num*D)