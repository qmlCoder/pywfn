"""
计算DFT格点
"""
from pywfn.base import Mol
from pywfn.data.elements import elements
from pywfn.data import lebedev
import numpy as np

class Calculator:
    def __init__(self,mol:Mol) -> None:
        self.mol=mol
        self.Nrad=[30,45,60]

    def radGrid(self,atmic:int):
        pi=np.pi
        nrad=self.Nrad[0]
        nrad=80
        R=elements[atmic].radius
        # print(f'{R=}')
        if atmic!=1:R/=2.
        rs=[]
        ws=[]
        for i in range(1,nrad+1):
            xi=np.cos(pi*i/(nrad+1))
            ri=R*(1.+xi)/(1.-xi)
            wi=2.*pi/(nrad+1)*R**3.*(1.+xi)**2.5/(1.-xi)**3.5
            # xi=i/(nrad+1)
            # ri=R*xi**2/(1.-xi)**2
            # wi=2.*R**3/(nrad+1)*xi**5/(1.-xi)**7
            # print(f'{i},{i/(nrad+1)},{xi:>10.4f},{ri:>10.4f},{wi:>10.4f}')
            if ri>15:continue
            rs.append(ri)
            ws.append(wi)
        # print('ri',len(rs))
        rs=np.array(rs,dtype=np.float32)
        ws=np.array(ws,dtype=np.float32)*4.*pi
        return rs,ws

    def sphGrid(self):
        result=lebedev.LD0074()
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
            
