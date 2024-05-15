"""
键级也不止一种，都在这里实现吧
"""
from pywfn.base import Mol
from pywfn.atomprop import direction
from pywfn.maths import CM2PM
from pywfn.utils import printer
import numpy as np
from itertools import product

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
        dirCaler=direction.Calculator(self.mol)
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
        dirCaler=direction.Calculator(self.mol)
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
        values=np.sqrt(result[:,-1])
        indexs=values>0
        result[:,-1]=values
        return result[indexs,:]

    def hmo(self):
        self.bond:list[int]=None
        # 1.建立系数矩阵
        atms=self.mol.heavyAtoms
        natm=len(atms)
        BM=np.zeros(shape=(natm,natm)) # 键连矩阵
        DM=np.zeros_like(BM) # 键长矩阵
        for i,j in product(range(natm),range(natm)):
            a1,a2=atms[i],atms[j]
            if a1>=a2:continue
            bond=self.mol.atom(a1).coord-self.mol.atom(a2).coord
            dist=np.linalg.norm(bond)
            DM[i,j]=dist
            DM[j,i]=dist
            if dist>1.7*1.889:continue
            BM[i,j]=1.0
            BM[j,i]=1.0
        # 2.求解
        e,C=np.linalg.eigh(BM) # 矩阵对角化
        nele=int(len(atms)-self.mol.charge) #电子数量
        idxs=np.argsort(e)[:nele//2] # 占据轨道
        CM=C[:,idxs] # 每一列对应一个特征向量
        # 3.构建键级矩阵
        result=[]
        OM=np.zeros_like(BM)
        for i,j in product(range(natm),range(natm)):
            order=np.sum(CM[i,:]*CM[j,:])*2
            OM[i,j]=order
            if i>=j:continue
            if DM[i,j]>1.7*1.889:continue
            result.append([atms[i],atms[j],order])
        result=np.array(result)
        # vals=result[:,-1]
        # vals=np.sqrt(vals**2)
        # result[:,-1]=vals
        return np.abs(result)

    def multiCenter(self,atms:list[int]):
        """计算多中心键级"""
        pass