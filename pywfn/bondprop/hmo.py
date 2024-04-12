"""
使用HMO方法计算键级
"""
from pywfn.bondprop import Caler
import numpy as np
from pywfn.base import Mol
from pywfn.utils import printer
class Calculator(Caler):
    def __init__(self,mol:"Mol"):
        self.mol=mol
        self.bond:list[int]=None
        # 1.建立系数矩阵
        atoms=self.mol.heavyAtoms
        self.atomMap={atom.idx:i for i,atom in enumerate(atoms)} #根据真实的原子序数，索引矩阵元素
        atomNum=len(atoms)
        self.M=np.zeros(shape=(atomNum,atomNum))
        for i,a in enumerate(atoms):
            for j,b in enumerate(atoms):
                dist=np.linalg.norm(a.coord-b.coord)
                if dist<1.7 and i!=j:
                    self.M[i,j]=1
        # 2.求解
        e,C=np.linalg.eigh(self.M) # 矩阵对角化
        eleNum=int(len(atoms)-self.mol.charge) #电子数量
        if (eleNum-2)%4!=0:
            printer.warn('该分子电子不满足4n+2规则')
        idxs=np.argsort(e)[:eleNum//2] # 占据轨道
        # 3.构建键级矩阵
        self.O=np.zeros_like(self.M)
        printer.info('列出所有键级绝对值大于0.1的结果：')
        for i,a in enumerate(atoms):
            for j,b in enumerate(atoms):
                order=np.sum(C[i,idxs]*C[j,idxs])*2
                self.O[i,j]=order
                dist=np.linalg.norm(a.coord-b.coord)
                if i<j and dist>1.7:
                    printer.log(f'{atoms[i].idx}-{atoms[j].idx}: {order:.4f}')
        pass
    
    def calculate(self):
        idx1,idx2=self.bond
        if idx1 not in self.atomMap.keys():return 0
        if idx2 not in self.atomMap.keys():return 0
        return self.O[idx1-1,idx2-1]