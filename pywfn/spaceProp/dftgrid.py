"""
计算DFT格点
"""
from pywfn.base import Mol
from pywfn.data.elements import elements
from pywfn.data import lebedev

import numpy as np
import sys

class Calculator:
    def __init__(self,mol:Mol) -> None:
        self.mol=mol
        self.nrad=74 # 径向格点数量
        self.nsph=86 # 球面格点数量

    def radGrid(self,atmic:int):
        pi=np.pi
        R=elements[atmic].radius
        if atmic!=1:R/=2.
        R=1.0
        rs=[]
        ws=[]
        icut=0
        for i in range(1,self.nrad+1):
            xi=np.cos(i*pi/(self.nrad+1))
            ri=R*(1. + xi)/(1. - xi)
            wi=2*pi/(self.nrad+1)*R**3*(1+xi)**2.5/(1.-xi)**3.5*4*pi
            # if ri>15:continue # 径向截断距离
            rs.append(ri)
            ws.append(wi)
            if ri<10 and icut==0:
                icut=i-1
            # print(f'radPos,{i:>10}{ri:>20.4f}{wi:>20.4f}')
        rs=np.array(rs)
        ws=np.array(ws)
        cuts=np.ones(self.nrad)
        cuts[:icut]=0.
        # idxs=np.argwhere(rs>10)
        # cuts[idxs]=0. # 根据multiwfn中的径向截断点得到的
        # ws[idxs]=0.   # 根据multiwfn中的径向截断点得到的
        # for i in range(self.nrad):
        #     print(f'radpos,{rs[i]:>20.4f}{ws[i]:>20.4f}{cuts[i]:>6.0f}')
        return rs,ws,cuts

    def sphGrid(self)->tuple[np.ndarray,np.ndarray]:
        """计算lebdev球面格点

        Returns:
            tuple[np.ndarray,np.ndarray]: 格点坐标,格点权重
        """
        match(self.nsph):
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
        # nsph=len(result)
        # for i in range(nsph):
        #     x,y,z,w=result[i]
        #     print(f'sphPos,{i+1:>5}{x:>10.4f}{y:>10.4f}{z:>10.4f}{w:>10.4f}')
        return result[:,:3],result[:,-1]

    def dftGrid(self,atm:int):
        """计算原子DFT格点"""
        cords=[]
        weits=[]
        gcuts=[]
        atom=self.mol.atom(atm)
        radRs,radWs,rcuts=self.radGrid(atom.atomic)
        sphCs,sphWs=self.sphGrid()
        # print('radRs',radRs.shape)
        # print('sphCs',sphCs.shape)
        # print(self.nrad,self.nsph)

        for rr,rw,ic in zip(radRs,radWs,rcuts):
            for sc,sw in zip(sphCs,sphWs):
                cords.append(sc*rr)
                weits.append(rw*sw)
                gcuts.append(ic)
        cords=np.array(cords,dtype=np.float64)
        weits=np.array(weits,dtype=np.float64)
        return cords,weits,gcuts
    
    def atmGrid(self,atm:int):
        """单个原子的网格点坐标，以原子为中心"""
        atmGrid,atmWeit,gcuts=self.dftGrid(atm)
        atmGrid=atmGrid+self.mol.atom(atm).coord.reshape(1,3) #移动到以原子为中心
        return atmGrid,atmWeit
    
    # def a2mGrid_(self,atm:int):
    #     """将原子中心坐标映射到分子中，改变权重"""
    #     from pywfn.maths import flib
    #     atmGrid,atmWeit=self.atmGrid(atm)
    #     nGrid=len(atmGrid)
    #     natm=len(self.mol.atoms)
    #     atmPos=self.mol.coords.copy() # 分子坐标
    #     atmRad=np.array(self.mol.atoms.radius)
    #     atmDis=self.mol.atoms.LM
    #     a2mGrid,a2mWeit=flib.a2mWeight(atm,nGrid,atmGrid,atmWeit,natm,atmPos,atmRad,atmDis)
    #     assert True not in np.isnan(atmWeit),"不应该有nan"
    #     return a2mGrid,a2mWeit
    
    def a2mGrid(self,atm:int):
        from pywfn.maths import rlib
        atmGrid,atmWeit=self.atmGrid(atm)
        atmPos=self.mol.coords.copy() # 分子坐标
        atmRad=self.mol.atoms.radius
        a2mWeits=rlib.a2m_weits_rs(atm-1,atmGrid.tolist(),atmWeit.tolist(),atmPos.tolist(),atmRad) # type: ignore
        a2mWeits=np.array(a2mWeits)
        return atmGrid,a2mWeits

    
    def molGrid(self):
        """整个分子的网格点坐标"""
        molGrid=[]
        molWeit=[]
        for atom in self.mol.atoms:
            a2mGrid,a2mWeit=self.a2mGrid(atom.idx)
            molGrid.append(a2mGrid)
            molWeit.append(a2mWeit)
        molPos=np.vstack(molGrid)
        molWei=np.concatenate(molWeit)
        return molPos,molWei