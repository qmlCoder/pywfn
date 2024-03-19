"""
计算原子p轨道的投影，作为新的指标
"""

from pywfn.base import Mol,Atom
import numpy as np

class Calculator:
    def __init__(self,mol:Mol) -> None:
        self.mol=mol
        self.vects:list[np.ndarray]=None # 选择的方向
        self.atoms:list[int]=None # 选择的原子

    def calculate(self):
        assert self.vects is not None,'未指定方向'
        assert self.atoms is not None,'未指定原子'
        # 获取原子的p轨道系数
        results=np.zeros(len(self.vects))
        for i,(idx,vect) in enumerate(zip(self.atoms,self.vects)):
            atom:Atom=self.mol.atom(idx)
            values=[]
            for obt in self.mol.O_obts:
                Cop=atom.get_pProj(vect,obt)
                assert len(Cop)!=0,'没有p轨道'
                values+=[np.linalg.norm(e) for e in Cop]
                # Cos=atom.get_sProj(self.vect,obt)
                # values+=[np.dot(e/np.sqrt(3)*self.vect,self.vect) for e in Cos]
            
            results[i]=np.mean(values)
        return results