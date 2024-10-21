from pywfn.base import Mol
from pywfn.spaceProp import density,dftgrid
from pywfn.maths import flib
import numpy as np


class Calculator:
    def __init__(self,mol:Mol) -> None:
        self.mol=mol
        self.cords=np.array([[0.,0.,0.]])
        self.grids,self.weits=dftgrid.Calculator(mol).molGrid()
    
    def set_grid(self,cords:np.ndarray): # 自己定义格点
        self.cords=cords

    def nucPotential_raw(self): # 计算原子势
        vals=np.zeros(shape=len(self.cords))
        for i,each in enumerate(self.cords):
            val=0.0
            for atom in self.mol.atoms:
                r=np.linalg.norm(each-atom.coord) # 两原子的距离
                if r<1e-6:continue
                val+=atom.atomic/r
            vals[i]=val
        return vals
    
    def nucPotential_lib(self)->np.ndarray: # 计算原子势
        nucs=np.array(self.mol.atoms.atomics,dtype=np.int64)
        vals=flib.nucPotential(self.cords,nucs,self.mol.coords)
        return vals
    
    # def elePotential_raw(self): # 计算电子势
    #     denCaler=density.Calculator(self.mol)
    #     denCaler.set_grid(self.grid)
    #     dens=denCaler.molDens_lib()
    #     res=np.zeros(len(self.grid))
    #     for i,r in enumerate(self.grid):
    #         val=0.0
    #         for j,r_ in enumerate(self.grid):
    #             if i==j:continue
    #             dist=np.linalg.norm(r-r_)
    #             if dist<1e-6:continue
    #             val+=dens[j]/dist
    #         res[i]=val
    #     return np.array(res)
    
    def elePotential_lib(self)->np.ndarray: # 计算电子势
        denCaler=density.Calculator(self.mol)
        denCaler.set_grid(self.grids)
        dens=denCaler.molDens_lib()
        vals=flib.elePotential(self.cords,self.grids,self.weits,dens)
        return vals
