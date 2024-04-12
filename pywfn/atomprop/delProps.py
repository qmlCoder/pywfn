"""
三个分子的差值确定的属性
- 福井函数
- parr函数
"""
from pywfn.base import Mol
from pywfn.atomprop import atomCharge,atomSpin,dirProps,AtomCaler
import numpy as np
from typing import Literal

class Calculator(AtomCaler):
    def __init__(self,mols:list[Mol]) -> None:
        self.mols=mols
        self.natm=len(mols[0].atoms)
        self.vects=None
        self.atoms=None
        self.chrg=None # 电荷类型
        self.prop=None # 计算属性
        self.func:Literal['fukui','parr']=None # 函数类型
        self.calers=[]
        
    
    def calculate(self)->np.ndarray:
        
        funcs={'fukui':self.fukuiFunc,'parr':self.parrFunc}
        if self.vects:
            assert self.vects is not None,'没有指定方向'
            self.calers=[dirProps.Calculator(mol) for mol in self.mols]
            for caler in self.calers:
                caler.atoms=self.atoms
                caler.vects=self.vects
                caler.chrg=self.chrg
                caler.prop=self.prop
            return funcs[self.func]
        else:
            if self.func=='fukui':
                self.calers=[atomCharge.Calculator(mol) for mol in self.mols]
            if self.func=='parr':
                self.calers=[atomSpin.Calculator(mol) for mol in self.mols]
            return funcs[self.func]

    
    def fukuiFunc(self)->np.ndarray:
        assert len(self.mols)==3,'需要三个分子'
        charges=np.zeros(shape=(self.natm,3))
        dcharge=np.zeros(shape=(self.natm,2))
        for c,caler in enumerate(self.calers):
            charge=caler.calculate()
            charges[:,c]=charge
        dcharge[:,0]=-(charges[:,0]-charges[:,1])
        dcharge[:,1]=-(charges[:,1]-charges[:,2])
        return dcharge
    
    def parrFunc(self)->np.ndarray:
        assert len(self.mols)==2,'需要两个分子'
        spins=np.zeros(shape=(self.natm,2))
        for c,caler in enumerate(self.calers):
            spin=caler.calculate()
            spins[:,c]=spin
        return spins
    