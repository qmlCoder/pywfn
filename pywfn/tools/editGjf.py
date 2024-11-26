"""
给gjf文件添加或减少一个电子
"""

from pywfn.base import Mol
from pywfn.writer import GjfWriter
from pathlib import Path
from pywfn.data.elements import elements
import numpy as np
from pywfn.shell import Shell
from pywfn.utils import printer

class Tool():
    def __init__(self,mol:Mol) -> None:
        self.mol=mol
        self.syms=list(mol.atoms.syms)
        self.xyzs=mol.atoms.xyzs
        self.nela=self.mol.eleNum[0]
        self.nelb=self.mol.eleNum[1]
        fold=Path(self.mol.reader.path).parent # 文件夹
        stem=Path(self.mol.reader.path).stem   # 文件名(无后缀)
        self.root=f'{fold}/{stem}'
    
    @property
    def multiply(self):
        atomic=sum([elements[sym].atomic for sym in self.syms]) #分子核电荷数
        charge=atomic-self.nela-self.nelb
        spin=self.nela-self.nelb+1
        return charge,spin
    
    def addElectron(self,delta:int):
        """分子加减电子，影响电荷，加一个电子电荷变为-1

        Args:
            delEle (int): 加减电子数
        """
        
        sat1=delta>0
        for i in range(abs(delta)):
            sat2=self.nela==self.nelb
            match (sat1,sat2):
                case (True,True):
                    self.nela+=1
                case (True,False):
                    self.nelb+=1
                case (False,True):
                    self.nelb-=1
                case (False,False):
                    self.nela-=1
    
    def addAtoms(self,syms:list[str],atms:np.ndarray):
        self.syms=self.syms+syms
        self.xyzs=np.concatenate([self.xyzs,atms],axis=0)
    
    def addRingBq(self,rings:list[list[int]]):
        xyzs=[]
        syms=['Bq']*len(rings)
        for ring in rings:
            xyz=self.xyzs[ring].copy().mean(axis=0)
            xyzs.append(xyz)
        self.addAtoms(syms,np.array(xyzs))
    
    def save(self,path:str):
        writer=GjfWriter()
        writer.charge,writer.spin=self.multiply
        writer.syms=self.syms
        writer.xyzs=self.xyzs
        writer.chk=Path(path).stem
        writer.save(path)