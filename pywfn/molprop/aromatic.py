"""
计算分子的芳香性
使用pi键级的标准差表示
"""
from pywfn.base import Mol
from pywfn.atomprop import charge
from pywfn.bondprop import order
import numpy as np

class Calculator:
    def __init__(self, mol:Mol):
        self.mol = mol
        self.caler=order.Calculator(mol)
        self.ratio=0.5
    
    def pisd_v1(self): # 版本1 直接用键级标准差
        result=self.caler.pi_pocv()
        orders=result[:,-1]
        forders=[] # 过滤掉C-H键
        for i in range(len(orders)):
            atm1=int(result[i,0])
            atm2=int(result[i,1])
            if self.mol.atom(atm1).atomic==1:continue
            if self.mol.atom(atm2).atomic==1:continue
            forders.append(orders[i])
        return np.std(forders).item()
    
    def pisd_v2(self): # 版本2 使用键级均值和标准差
        result=self.caler.pi_pocv()
        orders=result[:,-1]
        forders=[] # 过滤掉C-H键
        for i in range(len(orders)):
            atm1=int(result[i,0])
            atm2=int(result[i,1])
            if self.mol.atom(atm1).atomic==1:continue
            if self.mol.atom(atm2).atomic==1:continue
            forders.append(orders[i])
        mean=np.mean(forders)
        stds=np.std(forders)
        ratio=self.ratio
        return (ratio*mean-(1-ratio)*stds).item()

    def pimed(self): # 版本2 使用键级均值和标准差
        D=0.2
        result=self.caler.pi_pocv()
        orders=result[:,-1]
        forders=[] # 过滤掉C-H键
        for i in range(len(orders)):
            atm1=int(result[i,0])
            atm2=int(result[i,1])
            if self.mol.atom(atm1).atomic==1:continue
            if self.mol.atom(atm2).atomic==1:continue
            forders.append((orders[i]-0.66)**2)
        # mean=np.mean(forders)
        # stds=np.std(forders)
        # ratio=self.ratio
        nbond=self.mol.bonds.num
        return (1-np.sum(forders)/(nbond*D)).item()
    
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
        # print(vals)
        return 1-val/(self.mol.bonds.num*D)
