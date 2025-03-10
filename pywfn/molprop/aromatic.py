"""
计算分子的芳香性
使用pi键级的标准差表示
"""
from pywfn.base import Mol
from pywfn.atomprop import charge
from pywfn.bondprop import order as orderProp
import numpy as np
from pywfn.shell import Shell
from pywfn.utils import printer

class Calculator:
    def __init__(self, mol:Mol):
        self.mol = mol
        self.ratio=0.5
    
    def pisd(self,ring:list[int]|None=None): # 版本1 直接用键级标准差
        caler=orderProp.Calculator(self.mol)
        result=caler.pi_pocv()
        orders=result[:,-1]
        forders=[] # 过滤掉C-H键
        for i,order in enumerate(orders):
            atm1=int(result[i,0])
            atm2=int(result[i,1])
            if self.mol.atom(atm1).atomic==1:continue
            if self.mol.atom(atm2).atomic==1:continue
            if len(self.mol.atom(atm1).neighbors)==4:continue
            if len(self.mol.atom(atm2).neighbors)==4:continue
            if ring is None:
                print(f'{atm1:>2}-{atm2:>2}:{order:>10.4f}')
                forders.append(order)
            elif (atm1 in ring and atm2 in ring):
                print(f'{atm1:>2}-{atm2:>2}:{order:>10.4f}')
                forders.append(order)
            
        if len(forders)==0:
            printer.warn('没有键')
            return 0
        return np.std(forders).item()
    
    def pimsd(self,ring:list[int]|None=None,ratio=0.5): # 版本2 使用键级均值和标准差
        caler=orderProp.Calculator(self.mol)
        result=caler.pi_pocv()
        orders=result[:,-1]
        forders=[] # 过滤掉C-H键
        for i,order in enumerate(orders):
            atm1=int(result[i,0])
            atm2=int(result[i,1])
            if self.mol.atom(atm1).atomic==1:continue
            if self.mol.atom(atm2).atomic==1:continue
            if ring is None:
                print(f'{atm1:>2}-{atm2:>2}:{order:>10.4f}')
                forders.append(order)
            elif (atm1 in ring and atm2 in ring):
                print(f'{atm1:>2}-{atm2:>2}:{order:>10.4f}')
                forders.append(order)
        mean=np.mean(forders)
        stds=np.std(forders)
        return (ratio*mean-(1-ratio)*stds).item()

    def pimed(self): # 使用键级类比于HOMED方法
        D=0.2
        caler=orderProp.Calculator(self.mol)
        result=caler.pi_pocv()
        orders=result[:,-1]
        forders=[] # 过滤掉C-H键
        for i in range(len(orders)):
            atm1=int(result[i,0])
            atm2=int(result[i,1])
            if self.mol.atom(atm1).atomic==1:continue
            if self.mol.atom(atm2).atomic==1:continue
            forders.append((orders[i]-0.66)**2)
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