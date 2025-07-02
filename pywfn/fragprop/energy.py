from pywfn.base import Mole
from pywfn.maths import rlib
from pywfn.gridprop import dftgrid,density

import numpy as np

class Calculator:
    def __init__(self,mole:Mole) -> None:
        self.mole=mole

    def EIEBA(self,fragA:list[int],fragB:list[int]):
        gridCaler=dftgrid.Calculator(self.mole)
        densCaler=density.Calculator(self.mole)
        # grids,weits=gridCaler.molGrid()
        grids=[]
        weits=[]
        frags=[[],[]]
        rhos=[]
        for f,frag in enumerate([fragA,fragB]):
            fragGrids,fragWeits=gridCaler.fragGrid(frag)
            print(f,fragGrids.mean(axis=0))
            grids.append(fragGrids)
            weits.append(fragWeits)
            dens=densCaler.molDens(fragGrids,0)[0] # 电子密度在原子空间内的分布 / 原子电子密度在分子空间的分布
            rhos.append(dens)
            for atm in frag:
                atom=self.mole.atom(atm)
                px,py,pz=atom.coord
                frags[f].append([atom.atomic,px.item(),py.item(),pz.item()])
                
        gridsA=grids[0]
        gridsB=grids[1]
        weitsA=weits[0]
        weitsB=weits[1]
        # print(frags[0])
        # print(frags[1])
        # print(fragB)
        rhosA=rhos[0]
        rhosB=rhos[1]
        # print(rhosA.shape)
        # print(rhosB.shape)
        print(np.sum(rhosA*weitsA))
        print(np.sum(rhosB*weitsB))
        eleEN,eleEE,eleNN=rlib.calc_EIEBA_rs(frags[0],frags[1],rhosA,rhosB,gridsA,gridsB,weitsA,weitsB)
        return eleEN,eleEE,eleNN
        