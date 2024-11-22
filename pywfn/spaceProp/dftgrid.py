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
        self.nsph=434 # 球面格点数量

    def radGrid(self,atmic:int):
        pi=np.pi
        R=elements[atmic].radius
        # print(f'{R=}')
        if atmic!=1:R/=2.
        rs=[]
        ws=[]
        for i in range(1,self.nrad+1):
            xi=np.cos(pi*i/(self.nrad+1))
            ri=R*(1.+xi)/(1.-xi)
            wi=2.*pi/(self.nrad+1)*R**3.*(1.+xi)**2.5/(1.-xi)**3.5
            if ri>15:continue
            rs.append(ri)
            ws.append(wi)
        rs=np.array(rs)
        ws=np.array(ws)*4.*pi
        return rs,ws

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
            case 434:result=lebedev.LD0434()
            case _:
                print("原子角度格点数量不匹配!!")
                sys.exit(1)
        return result[:,:3],result[:,-1]

    def dftGrid(self,atm:int):
        """计算原子DFT格点"""
        cords=[]
        weits=[]
        atom=self.mol.atom(atm)
        radRs,radWs=self.radGrid(atom.atomic)
        sphCs,sphWs=self.sphGrid()
        for rr,rw in zip(radRs,radWs):
            for sc,sw in zip(sphCs,sphWs):
                cords.append(sc*rr)
                weits.append(rw*sw)
        cords=np.array(cords,dtype=np.float32)
        weits=np.array(weits,dtype=np.float32)
        return cords,weits
    
    def atmGrid(self,atm:int):
        """单个原子的网格点坐标，以原子为中心"""
        atmGrid,atmWeit=self.dftGrid(atm)
        atmGrid=atmGrid+self.mol.atom(atm).coord.reshape(1,3) #移动到以原子为中心
        return atmGrid,atmWeit
    
    def a2mGrid(self,atm:int):
        """将原子中心坐标映射到分子中，改变权重"""
        from pywfn.maths import flib
        atmGrid,atmWeit=self.atmGrid(atm)
        nGrid=len(atmGrid)
        natm=len(self.mol.atoms)
        atmPos=self.mol.coords.copy() # 分子坐标
        atmRad=np.array(self.mol.atoms.radius)
        atmDis=self.mol.atoms.LM
        a2mGrid,a2mWeit=flib.a2mWeight(atm,nGrid,atmGrid,atmWeit,natm,atmPos,atmRad,atmDis)
        assert True not in np.isnan(atmWeit),"不应该有nan"
        return a2mGrid,a2mWeit
    
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