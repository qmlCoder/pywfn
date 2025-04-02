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
from pywfn.data import consts

class Calculator:
    def __init__(self, mol:Mole):
        self.mol = mol
        self.ratio=0.5
    
    def pisd(self,ring:list[int]|None=None): # 版本1 直接用键级标准差
        caler=orderProp.Calculator(self.mol)
        bonds,orders=caler.piOrder_pocv()
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
    
    def pimsd(self,ring:list[int]|None=None,ratio=0.5) -> float: # 版本2 使用键级均值和标准差
        caler=orderProp.Calculator(self.mol)
        bonds,orders=caler.piOrder_pocv()
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

    def pimed(self) -> float: # 使用键级类比于HOMED方法
        D=0.2
        caler=orderProp.Calculator(self.mol)
        bonds,orders=caler.piOrder_pocv()
        forders=[] # 过滤掉C-H键
        for i,((a1,a2),order) in enumerate(zip(bonds,orders)):
            if self.mol.atom(a1).atomic==1:continue
            if self.mol.atom(a2).atomic==1:continue
            forders.append((order-0.66)**2)
        nbond=self.mol.bonds.num
        result=1-sum(forders)/(nbond*D)
        return result
    
    def homed(self) -> float:
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
    
    def HOMER(self,rings:list[list[int]]):
        """
        Phys. Chem. Chem. Phys., 2023, 25, 16763
        """
        paras={
            'CC':[1.437,950.74],
            'CN':[1.390,506.43],
            'NC':[1.390,506.43],
            'NN':[1.375,187.36],
            'CO':[1.379,164.96],
            'OC':[1.379,164.96],
        }
        vals:list[float]=[]
        for ring in rings:
            val=0.0
            n=0
            for bond in self.mol.bonds:
                a1,a2=bond.ats
                if a1 not in ring or a2 not in ring:continue # 两个原子都必须在环中
                s1=self.mol.atom(a1).sym
                s2=self.mol.atom(a2).sym
                key=f'{s1}{s2}'
                assert key in paras.keys(),"只能包含C,N,O三种元素"
                # print(bond,bond.length*consts.Bohr)
                n+=1
                Ropt,alp=paras[key]
                val+=alp*(bond.length*consts.Bohr-Ropt)**2
            # print(n)
            vals.append(1-val/n)
        return vals


