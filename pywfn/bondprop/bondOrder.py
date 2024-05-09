"""
键级也不止一种，都在这里实现吧
"""
from pywfn.base import Mol
from pywfn.atomprop import atomDirect

import numpy as np

class Calculator:
    def __init__(self,mol:Mol) -> None:
        self.mol=mol

    def mayer(self)->np.ndarray:
        """计算mayer键级"""
        # 获取密度矩阵 P
        PM=self.mol.PM
        # 获取重叠矩阵
        SM=self.mol.SM
        PS=PM@SM
        OM=PS*PS.T
        orders=[]
        for bond in self.mol.bonds:
            u1,l1=bond.a1.obtBorder
            u2,l2=bond.a2.obtBorder
            order=np.sum(OM[u1:l1,u2:l2])
            orders.append([bond.atm1,bond.atm2,order])
        return np.array(orders)
    
    def dirMayer(self,bonds:list[list[int]]):
        """带有方向的Mayer键级[d,6](a1,a2,x,y,z,v)"""
        dirCaler=atomDirect.Calculator(self.mol)
        keys=self.mol.bonds.keys
        obts=self.mol.O_obts
        CMo=self.mol.CM
        result=[]
        for a1,a2 in bonds:
            dirs=dirCaler.reaction(a1)
            atms=[a1]*len(dirs)
            if a1>a2:a1,a2=a2,a1
            key=f'{a1}-{a2}'
            idx=keys.index(key)
            for d,(atm,dir_) in enumerate(zip(atms,dirs)):
                CMp=self.mol.projCM(obts,[atm],[dir_],False,True,False)
                self.mol.props['CM']=CMp
                a1_,a2_,order=self.mayer()[idx]
                assert a1==a1_ and a2==a2_,"原子不对应"
                x,y,z=dir_
                result.append([a1,a2,x,y,z,order])
        self.mol.props['CM']=CMo
        return np.array(result)
    
    def piOrder(self):
        """计算pi键级"""
        pass

    
    def multiCenter(self,atms:list[int]):
        """计算多中心键级"""
        pass