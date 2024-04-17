from pywfn.base import Mol,Atom
import numpy as np
from pywfn.data import Elements
from functools import lru_cache
elements=Elements()
from pywfn.utils import printer
from pywfn.atomprop import lutils,AtomCaler
from pywfn import maths
from typing import Literal

class Calculator(AtomCaler):
    def __init__(self,mol:"Mol"):
        self.logTip:str=''
        self.mol=mol
        self.chrg:Literal['mulliken','lowdin']='mulliken'
    
    def calculate(self)->np.ndarray:
        if self.chrg=='mulliken':
            return self.mulliken()
        if self.chrg=='lowdin':
            return self.lowdin()
    
    def mulliken(self):
        # 计算密度矩阵
        self.logTip='Mulliken电荷分布'
        PM=self.mol.PM
        # PS=PM*mol.SM
        # PSS=PS.sum(axis=0)
        # 矩阵乘法的迹的加和=矩阵对应元素乘积之和
        PS=PM@self.mol.SM
        EV=np.diagonal(PS)
        atoms=self.mol.atoms
        charges=np.zeros(len(atoms))
        for a,atom in enumerate(atoms):
            a1,a2=atom.obtBorder
            elect=EV[a1:a2].sum()
            charges[a]=atom.atomic-elect
        return charges

    def lowdin(self):
        """
        计算每个原子的lowdin电荷
        """
        self.logTip='Lowdin电荷分布'
        PM=self.mol.PM
        SM=self.mol.SM
        # 计算矩阵的1/2次方
        v, Q=np.linalg.eig(SM)
        V=np.diag(v)
        V_=V**0.5
        Q_=np.linalg.inv(Q)
        SM_half=Q@(V_@Q_)
        SPS=SM_half@(PM@SM_half)
        eleNums=np.diagonal(SPS)
        charges=np.zeros(len(self.mol.atoms))
        for a,atom in enumerate(self.mol.atoms):
            u,l=atom.obtBorder
            eleNum=eleNums[u:l].sum()
            charges[a]=atom.atomic-eleNum
        return charges
    
    def hirshfeld(self):
        """
        计算原子的Hirshfeld电荷
        """
        self.mol.bohr=True
        from pywfn.data import sphGrid,radDens
        coords=sphGrid.gridData[:,:3]# 原点为0的坐标
        weight=sphGrid.gridData[:,3]
        obts=self.mol.O_obts
        atoms:list[int]=self.mol.atoms.indexs
        chargs=np.zeros(shape=(len(atoms)))
        for a1,atom1 in enumerate(self.mol.atoms): # 计算每一个原子的电荷
            coord=atom1.coord+coords # 该原子周围点的空间坐标
            # molDens=self.mol.get_dens(atoms,obts,coord)*weight # 计算分子的电子密度
            molDens=np.zeros(shape=(len(coords))) # 计算分子的电子密度
            proDens=np.zeros(shape=(len(coords)))
            for a2,atom2 in enumerate(self.mol.atoms): #计算前置分子密度
                radius=np.linalg.norm(coord-atom2.coord,axis=1) #所有坐标对应的半径
                dens1=radDens.get_radDens(atom2.atomic,radius)*weight # 插值法计算提前存储好的密度
                # print(a1==a2,a2+1,np.sum(dens))
                proDens+=dens1
                if a2==a1:atmDens=dens1
                dens2=atom2.get_dens(obts,coord-atom2.coord,weight)
                molDens+=dens2
                print(f'{dens1.sum():.4f},{dens2.sum():.4f}')
            ratio=np.divide(atmDens,proDens,out=np.zeros_like(atmDens),where=proDens!=0)
            chargs[a1]=np.sum(ratio*molDens)
            atmQ,proQ,molQ=np.sum(atmDens),np.sum(proDens),np.sum(molDens)
            print(f'{atmQ=:.4f},{proQ=:.4f},{molQ=:.4f}')
        return chargs
    
    def resStr(self)->str:
        """获取结果的打印内容"""
        satoms=lutils.atomIdxs(self.mol.atoms)
        charges=self.calculate()
        return lutils.atomValueStr(self.mol,satoms,charges)