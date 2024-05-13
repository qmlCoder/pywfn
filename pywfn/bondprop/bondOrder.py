"""
键级也不止一种，都在这里实现吧
"""
from pywfn.base import Mol
from pywfn.atomprop import atomDirect
from pywfn.maths import CM2PM

import numpy as np

class Calculator:
    def __init__(self,mol:Mol) -> None:
        self.mol=mol

    def mayer(self,*,PM=None,bonds=None)->np.ndarray:
        """计算指定键的mayer键级，可以有多个"""
        # 获取密度矩阵 P
        if PM is None:
            PM=self.mol.PM
        # 获取重叠矩阵
        SM=self.mol.SM
        PS=PM@SM
        OM=PS*PS.T
        orders=[]
        if bonds is None:
            bonds=[bond.ats for bond in self.mol.bonds]
        for a1,a2 in bonds:
            u1,l1=self.mol.atom(a1).obtBorder
            u2,l2=self.mol.atom(a2).obtBorder
            order=np.sum(OM[u1:l1,u2:l2])
            orders.append([a1,a2,order])
        order = np.array(orders)
        # print(order)
        return order
    
    def dirMayer(self,bonds:list[list[int,int]])->np.ndarray:
        """带有方向的Mayer键级[d,6](a1,a2,x,y,z,v)"""
        dirCaler=atomDirect.Calculator(self.mol)
        obts=self.mol.O_obts
        result=[]
        for a1,a2 in bonds:
            dirs=dirCaler.reaction(a1)
            atms=[a1]*len(dirs)
            if a1>a2:a1,a2=a2,a1
            for d,(atm,dir_) in enumerate(zip(atms,dirs)):
                CMp=self.mol.projCM(obts,[atm],[dir_],True,True)
                PMp=CM2PM(CMp,self.mol.O_obts,self.mol.oE)
                a1_,a2_,order=self.mayer(PM=PMp,bonds=[[a1,a2]]).flatten()
                assert a1==a1_ and a2==a2_,"原子不对应"
                x,y,z=dir_
                result.append([a1,a2,x,y,z,order])
        return np.array(result)
    
    def piOrder(self):
        """
        计算pi键级，每一个可能的π键计算出一个π键级
        """
        dirCaler=atomDirect.Calculator(self.mol)
        atms=[]
        dirs=[]
        for atom in self.mol.atoms:
            normal=dirCaler.normal(atom.idx) # 原子的法向量
            if normal is None:continue
            atms.append(atom.idx)
            dirs.append(normal)
        PMp=self.mol.projCM(self.mol.O_obts,atms,dirs,False,False)
        PMp=CM2PM(PMp,self.mol.O_obts,self.mol.oE)
        result=self.mayer(PM=PMp)
        result[:,-1]=np.sqrt(result[:,-1])
        return result

    
    def multiCenter(self,atms:list[int]):
        """计算多中心键级"""
        pass