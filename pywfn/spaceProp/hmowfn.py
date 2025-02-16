"""
休克尔分子轨道法轨道的可视化

将空间格点转换为原子局部坐标系的点，计算这些点上的波函数值
"""

from pywfn.base import Mol
from pywfn.reader import LogReader
from pywfn.atomprop import direction
from pywfn.maths import vector_angle
from pywfn.maths.mol import hmo

import numpy as np

def hmoWfn(Z:int,grid:np.ndarray):
    a0=0.589
    n=1/(4*np.sqrt(2*np.pi))*(Z/a0)**(5/2)
    r=np.linalg.norm(grid,axis=1) # 到原子核的距离
    wfn=n*np.exp(-Z*r*0.5/a0)*grid[:,-1]
    return wfn


class Calculator:
    def __init__(self,mol:Mol):
        self.mol=mol
        BM,es,CM,occs=hmo(self.mol)
        self.CM=CM

    def ato_wfns(self,grids:np.ndarray):
        """计算所有的原子轨道波函数
        """
        from pywfn.atomprop import direction
        dirCaler=direction.Calculator(self.mol)
        base=dirCaler.hmoBases()
        atms=self.mol.heavyAtoms
        natm=len(atms)
        wfns=np.zeros(shape=(natm,len(grids)))
        for i,atm in enumerate(atms):
            atom=self.mol.atom(atm)
            grid=(grids-atom.coord)@base[atm] # 将空间格点转为以原子为中心的坐标
            wfn=hmoWfn(atom.atomic,grid)
            wfns[i]=wfn
        return wfns
    
    def ato_wfn(self,grids,atm:int):
        from pywfn.atomprop import direction
        dirCaler=direction.Calculator(self.mol)
        base=dirCaler.hmoBases()
        atom=self.mol.atom(atm)
        grid=(grids-atom.coord)@base[atm] # 将空间格点转为以原子为中心的坐标
        wfn=hmoWfn(atom.atomic,grid)
        return wfn

    
    def obtWfn(self,grids:np.ndarray,obt:int)->np.ndarray:
        wfns=self.ato_wfns(grids)
        atms=self.mol.heavyAtoms
        natm=len(atms)
        wfn=np.zeros(len(grids))
        for i,atm in enumerate(atms):
            coe=self.CM[i,obt]
            wfn+=coe*wfns[i]
        return wfn