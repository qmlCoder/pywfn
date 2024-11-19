from pywfn.base import Mol
from pywfn.spaceProp import density,dftgrid
from pywfn.maths import flib
import numpy as np


class Calculator:
    def __init__(self,mol:Mol) -> None:
        self.mol=mol
        gridCaler=dftgrid.Calculator(mol)
        gridCaler.nrad=50
        gridCaler.nsph=74
        grids,weits=gridCaler.molGrid()
        self.grids=grids
        self.weits=weits

    def nucPotential_raw(self,grid:np.ndarray): # 计算原子势
        vals=np.zeros(shape=len(grid))
        for i,each in enumerate(grid):
            val=0.0
            for atom in self.mol.atoms:
                r=np.linalg.norm(each-atom.coord) # 两原子的距离
                if r<1e-6:continue
                val+=atom.atomic/r
            vals[i]=val
        return vals
    
    def nucPotential_lib(self,grid:np.ndarray)->np.ndarray: # 计算原子势
        nucs=np.array(self.mol.atoms.atomics,dtype=np.int64)
        vals=flib.nucPotential(grid,nucs,self.mol.coords)
        return vals
    
    def elePotential_lib(self,grid:np.ndarray)->np.ndarray: # 计算电子势
        denCaler=density.Calculator(self.mol)
        dens=denCaler.molDens_lib(self.grids) # 电子密度
        vals=flib.elePotential(grid,self.grids,self.weits,dens)
        return vals
