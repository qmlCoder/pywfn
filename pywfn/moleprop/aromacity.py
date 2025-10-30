"""
计算分子的芳香性
使用pi键级的标准差表示
"""
from pywfn.base.mole import Mole
from pywfn.atomprop import charge
from pywfn.bondprop import order as orderProp
import numpy as np
from pywfn import core
from pywfn.utils import printer
from pywfn.data import consts

class Calculator:
    def __init__(self, mole:Mole):
        self.mole = mole
        self.ratio=0.5
        self.caler=core.moleprop.aromacity.Calculator(mole.mole)  # type: ignore
    
    def PISD(self,rings:list[list[int]]|None=None): # 版本1 直接用键级标准差
        print(self.caler.PISD.__doc__)
        return self.caler.PISD(rings)
    
    # def PIMSD(self,ring:list[int]|None=None,ratio=0.5) -> float: # 版本2 使用键级均值和标准差
    #     caler=orderProp.Calculator(self.mole)
    #     bonds,orders=caler.pi_pocv()
    #     forders=[] # 过滤掉C-H键
    #     for i,((a1,a2),order) in enumerate(zip(bonds,orders)):
    #         if self.mole.atom(a1).atomic==1:continue
    #         if self.mole.atom(a2).atomic==1:continue
    #         if ring is None:
    #             print(f'{a1:>2}-{a2:>2}:{order:>10.4f}')
    #             forders.append(order)
    #         elif (a1 in ring and a2 in ring):
    #             print(f'{a1:>2}-{a2:>2}:{order:>10.4f}')
    #             forders.append(order)
    #     mean=np.mean(forders)
    #     stds=np.std(forders)
    #     return (ratio*mean-(1-ratio)*stds).item()

    # def PIMED(self) -> float: # 使用键级类比于HOMED方法
    #     D=0.2
    #     caler=orderProp.Calculator(self.mole)
    #     bonds,orders=caler.pi_pocv()
    #     forders=[] # 过滤掉C-H键
    #     for i,((a1,a2),order) in enumerate(zip(bonds,orders)):
    #         if self.mole.atom(a1).atomic==1:continue
    #         if self.mole.atom(a2).atomic==1:continue
    #         forders.append((order-0.66)**2)
    #     nbond=self.mole.bonds.num
    #     result=1-sum(forders)/(nbond*D)
    #     return result
    
    def HOMED(self) -> float:
        D=0.2
        idea=1.39645
        vals=[]
        for bond in self.mole.bonds:
            atm1,atm2=bond.atm1,bond.atm2
            if self.mole.atom(atm1).atomic==1:continue
            if self.mole.atom(atm2).atomic==1:continue
            vals.append((bond.length/1.889-idea)**2)
        val=sum(vals)
        return 1-val/(self.mole.bonds.num*D)

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
            for bond in self.mole.bonds:
                a1,a2=bond.ats
                if a1 not in ring or a2 not in ring:continue # 两个原子都必须在环中
                s1=self.mole.atom(a1).sym
                s2=self.mole.atom(a2).sym
                key=f'{s1}{s2}'
                assert key in paras.keys(),"只能包含C,N,O三种元素"
                # print(bond,bond.length*consts.Bohr)
                n+=1
                Ropt,alp=paras[key]
                val+=alp*(bond.length*consts.Bohr-Ropt)**2
            # print(n)
            vals.append(1-val/n)
        return vals


