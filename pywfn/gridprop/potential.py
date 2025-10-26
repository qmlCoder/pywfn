from pywfn.base.mole import Mole
from pywfn.gridprop import density,dftgrid
from pywfn import gridprop
import numpy as np
import json
from pathlib import Path


class Calculator(gridprop.SpaceCaler):
    def __init__(self,mole:Mole) -> None:
        self.mole=mole
        gridCaler=dftgrid.Calculator(mole)
        gridCaler.nrad=50
        gridCaler.fsph=74
        grids,weits=gridCaler.mol_grids()
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
    

    def nucPotential(self,qpos:np.ndarray):
        from pywfn import core
        nucs=np.array(self.mole.atoms.atomics,dtype=np.float64)
        xyzs=self.mole.xyzs
        vals=core.space.nuc_potential(qpos,xyzs,nucs) # type: ignore
        vals=np.array(vals)
        return vals
    
    def elePotential(self,qpos:np.ndarray):
        from pywfn import core
        densCaler=density.Calculator(self.mol)
        dens,_,_=densCaler.molDens(self.grids,level=0) # 电子密度
        vals=core.space.ele_potential(qpos,self.grids,self.weits,dens) # type: ignore
        vals=np.array(vals)
        return vals
    
    def molPotential(self,grid:np.ndarray)->np.ndarray: # 计算分子势
        nucPot=self.nucPotential(grid)
        elePot=self.elePotential(grid)
        return nucPot+elePot
    
    def electSurface(self)->tuple[np.ndarray,np.ndarray,np.ndarray]: # 计算电子势能面，网格是根据电子密度等值面确定的
        """计算电子势能面，网格是根据电子密度等值面确定的

        Returns:
            tuple[np.ndarray,np.ndarray,np.ndarray]: 顶点，面，静电势
        """
        densCaler=density.Calculator(self.mol)
        verts,faces=densCaler.vandSurf()
        elec=self.elePotential(verts)
        return verts,faces,elec
