"""
三个分子的差值确定的属性
- 福井函数
- parr函数
"""
from pywfn.base import Mol
from pywfn.atomprop import atomCharge,atomSpin,dirProps,AtomCaler,atomEnergy
import numpy as np
from typing import Literal

class Calculator(AtomCaler):
    def __init__(self,mols:list[Mol]) -> None:
        self.mols=mols
        self.natm=len(mols[0].atoms)
        self.vects=None # 可以计算方向属性
        self.atoms=None
        self.chrg=None # 电荷类型
        self.prop=None # 计算属性
        self.func:Literal['fukui','parr']=None # 函数类型
        self.calers=[]
        self.hists={}
        
    
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

    def setProject(self,stype:bool):
        """
        检查是否需要投影
        stype,修改(T)或撤回修改(F)
        """
        if self.vects is None:return
        for m,mol in enumerate(self.mols):
            if stype:
                CMo=mol.CM
                self.hists[f'CM_{m}']=CMo # 将老的存储下来
                CMn=mol.projCM(self.atoms,mol.O_obts,self.vects,True,False,False)
                mol.datas['CM']=CMn
            else:
                CMo=self.hists[f'CM_{m}'] # 恢复成老的
                mol.datas['CM']=CMo
    
    def fukuiFunc(self)->np.ndarray:
        self.setProject(True)
        assert len(self.mols)==3,'需要三个分子'
        charges=np.zeros(shape=(self.natm,3))
        dcharge=np.zeros(shape=(self.natm,2))
        for c,caler in enumerate(self.calers):
            charge=caler.calculate()
            charges[:,c]=charge
        dcharge[:,0]=-(charges[:,0]-charges[:,1])
        dcharge[:,1]=-(charges[:,1]-charges[:,2])
        self.setProject(False)
        return dcharge
    
    def parrFunc(self)->np.ndarray:
        self.setProject(True)
        assert len(self.mols)==2,'需要两个分子'
        spins=np.zeros(shape=(self.natm,2))
        for c,caler in enumerate(self.calers):
            spin=caler.calculate()
            spins[:,c]=spin
        self.setProject(False)
        return spins
    
    def atomEngs(self)->np.ndarray:
        """计算原子轨道能"""
        self.setProject(True)
        self.calers=[atomEnergy.Calculator(mol) for mol in self.mols]
        atmNum=len(self.mols[0].atoms)
        allEngs=np.zeros(shape=(atmNum,3))
        for c,caler in enumerate(self.calers):
            engs=caler.calculate()
            allEngs[:,c]=engs
        self.allEngs=allEngs
        delEngs=np.zeros(shape=(atmNum,2))
        delEngs[:,0]=allEngs[:,0]-allEngs[:,1]
        delEngs[:,1]=allEngs[:,1]-allEngs[:,2]
        self.setProject(False)
        return delEngs
    