"""
计算DFT格点
"""
from pywfn.base.mole import Mole
from pywfn.data.elements import elements
from pywfn.data import lebedev

import numpy as np
import sys
from pywfn import core

class Calculator:
    def __init__(self,mole:Mole) -> None:
        self.mole=mole
        self.nrad=99 # 径向格点数量
        self.fsph=74 # 球面格点数量
        self.caler=core.gridprop.dftgrid.Calculator(mole.mole) # type: ignore

    def rad_grids(self,atmic:int):
        R=elements[atmic].radius
        return self.caler.rad_grids(R,self.nrad)

    def sph_grids(self)->tuple[np.ndarray,np.ndarray]:
        """计算lebdev球面格点

        Returns:
            tuple[np.ndarray,np.ndarray]: 格点坐标,格点权重
        """
        match(self.fsph):
            case 6: result=lebedev.LD0006()
            case 14:result=lebedev.LD0014()
            case 26:result=lebedev.LD0026()
            case 74:result=lebedev.LD0074()
            case 86:result=lebedev.LD0086()
            case 230:result=lebedev.LD0230()
            case 434:result=lebedev.LD0434()
            case _:
                print("原子角度格点数量不匹配!!")
                sys.exit(1)
        return result[:,:3],result[:,-1]
    
    def atm_grids(self,iatm:int):
        return self.caler.atm_grids(iatm,self.nrad,self.fsph)
    
    def mol_grids(self):
        """整个分子的网格点坐标"""
        grids,weits=self.caler.mol_grids(self.nrad,self.fsph)
        grids=np.array(grids)
        weits=np.array(weits)
        return grids,weits
    
    def fragGrid(self,frag:list[int]):
        """分子中某个片段的格点"""
        return self.caler.frag_grids(frag)

    # def dftGrid(self,atm:int):
    #     """计算原子DFT格点"""
    #     cords=[]
    #     weits=[]
    #     gcuts=[]
    #     atom=self.mole.atom(atm)
    #     radRs,radWs,rcuts=self.radGrid(atom.atomic)
    #     sphCs,sphWs=self.sphGrid()

    #     for rr,rw,ic in zip(radRs,radWs,rcuts):
    #         for sc,sw in zip(sphCs,sphWs):
    #             cords.append(sc*rr)
    #             weits.append(rw*sw)
    #             gcuts.append(ic)
    #     cords=np.array(cords,dtype=np.float64)
    #     weits=np.array(weits,dtype=np.float64)
    #     return cords,weits,gcuts
    
    
    # def a2mGrid(self,atm:int):
    #     from pywfn import core
    #     atmGrid,atmWeit=self.atmGrid(atm)
    #     atmPos=self.mole.xyzs.copy() # 分子坐标
    #     atmRad=self.mole.atoms.radius
    #     a2mWeits=core.space.a2m_weits(atm-1,atmGrid.tolist(),atmWeit.tolist(),atmPos.tolist(),atmRad) # type: ignore
    #     a2mWeits=np.array(a2mWeits)
    #     return atmGrid,a2mWeits

    
