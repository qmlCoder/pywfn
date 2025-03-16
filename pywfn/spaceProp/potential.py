from pywfn.base import Mol
from pywfn.spaceprop import density,dftgrid
from pywfn import spaceprop
from pywfn.maths import flib,rlib
import numpy as np
import json
from pathlib import Path


class Calculator(spaceprop.SpaceCaler):
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
    
    def nucPotential_(self,grid:np.ndarray)->np.ndarray: # 计算原子势
        nucs=np.array(self.mol.atoms.atomics,dtype=np.int64)
        vals=flib.nucPotential(grid,nucs,self.mol.coords)
        return vals
    
    def nucPotential(self,qpos:np.ndarray):
        nucs=np.array(self.mol.atoms.atomics,dtype=np.float64)
        xyzs=self.mol.coords
        vals=rlib.nuc_potential_rs(qpos,xyzs,nucs) # type: ignore
        vals=np.array(vals)
        return vals
    
    # 静电式是根据电子密度计算出来的，所以控制不同原子的电子密度即可对应不同原子的静电式
    # 但是貌似不需要计算不同原子的静电式吧？需要的时候再说吧
    def elePotential_(self,grid:np.ndarray)->np.ndarray: # 计算电子势
        densCaler=density.Calculator(self.mol)
        dens,_,_=densCaler.molDens(self.grids,level=0) # 电子密度
        vals=flib.elePotential(grid,self.grids,self.weits,dens)
        return vals
    
    def elePotential(self,qpos:np.ndarray):
        densCaler=density.Calculator(self.mol)
        dens,_,_=densCaler.molDens(self.grids,level=0) # 电子密度
        vals=rlib.ele_potential_rs(qpos,self.grids,self.weits,dens) # type: ignore
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
