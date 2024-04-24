"""
轨道投影法+Mayer计算π键级
根据投影计算新的分子轨道系数矩阵,然后重构密度矩阵
"""

from pywfn.base import Mol
from typing import *
import numpy as np

from pywfn import maths
from pywfn.bondprop import Caler
from pywfn.utils import printer

class Calculator(Caler):
    def __init__(self,mol:Mol) -> None:
        self.mol=mol
        self.values={}
        self.vect:Union[None,np.ndarray]=None
        self.bond:list[int,int]=None
        self.zero:bool=True
        self.keep:bool=True
        self.ins:bool=True

    def calculate(self)->np.ndarray:
        """
        计算pi键级，当指定方向时计算方向键级
        """
        assert self.bond is not None,'没有指定键级'
        idx1,idx2=self.bond
        printer.info(f'计算{idx1} → {idx2}的键级')
        
        centerAtom=self.mol.atom(idx1)
        aroundAtom=self.mol.atom(idx2)
        if centerAtom.symbol=='H' or aroundAtom.symbol=='H':
            return 0.,0.
        obts=self.mol.O_obts #占据轨道的索引

        vect=self.vect
        if vect is None:
            normal=centerAtom.get_Normal(aroundAtom.idx)
        else:
            normal=vect/np.linalg.norm(vect)
        
        
        if len(centerAtom.neighbors)==3:
            printer.vector('法向方向: ',normal)
            return self.calerWay(obts,normal),0.
        elif len(centerAtom.neighbors) in [2,1]:
            bondVector=aroundAtom.coord-centerAtom.coord
            if normal is None:
                idxn=maths.search_sp2(idx1,self.mol)
                if idxn is None:
                    normal=centerAtom.get_vertObt(idx1,idx2)
                    printer.vector('轨道方向: ',normal)
                else:
                    normal=self.mol.atom(idxn).get_Normal(idx2)
                    printer.vector('相邻法向: ',normal)
            else:
                printer.vector('法向方向: ',normal)
            
            
            cross=np.cross(bondVector,normal)
            cross/=np.linalg.norm(cross)
            printer.vector('交叉方向: ',cross)
            vn,vc=self.calerWay(obts,normal),self.calerWay(obts,cross)
            return np.array([vn,vc])
        else:
            return np.array([0.,0.])

    def calerWay(self,obts:list[int],norm:np.ndarray):
        """计算某个方向的键级"""
        norm.flags.writeable=False
        length=np.linalg.norm(norm)
        assert abs(length-1)<1e-4,"向量长度应为1"
        idx1,idx2=self.bond
        """计算两个原子之间某个方向的键级"""
        atom1=self.mol.atom(idx1)
        atom2=self.mol.atom(idx2)
        a_1,a_2=atom1.obtBorder
        b_1,b_2=atom2.obtBorder
        atoms=[idx1,idx2]
        vects=[norm]*len(atoms)
        CM_=self.mol.projCM(atoms,obts,vects,zero=self.zero,keep=self.keep,ins=self.ins)

        oe=self.mol.oE
        SM=self.mol.SM
        
        PM_=maths.CM2PM(CM_,obts,oe)
        PS=PM_@SM
        order=np.sum(PS[a_1:a_2,b_1:b_2]*PS[b_1:b_2,a_1:a_2].T)
        return np.sqrt(np.abs(order))
    

    def resStr(self) -> str:
        order1,order2=self.calculate()
        return f'{order1:.4f},{order2:.4f}'
    
    def print(self):
        printer.res(self.resStr())