"""
各种带方向的原子属性放到一块吧
"""
from pywfn.base import Mol
from pywfn.atomprop import atomCharge, atomSpin,freeValence
from pywfn.atomprop import atomEnergy
from pywfn.atomprop import lutils
from pywfn.utils import printer
import numpy as np
from typing import Literal


class Calculator:
    """
    带有方向的各种原子属性的计算
    需要指定：电荷(chrg),prop(性质),atoms(原子),vects(方向)
    """
    def __init__(self,mol:Mol) -> None:
        self.mol=mol
        self.chrg:Literal['mulliken','lowdin']=None
        self.prop:Literal['charge','spin']=None
        self.vects:list[np.ndarray]=None
        self.atoms:list[int]=None
        self.zero:bool=True
        self.keep:bool=False
        self.abs:bool=False
        self.ins:bool=False

    def calculate(self):
        assert self.chrg is not None,'没有指定方法'
        assert self.prop is not None,'没有指定属性'
        assert self.vects is not None,'没有指定方向'
        assert self.atoms is not None,'没有指定原子'
        assert len(self.atoms)==len(self.vects),'方向与原子数量应该一致'
        # 获取投影后的系数轨道
        CM=self.mol.projCM(self.atoms,self.mol.O_obts,self.vects,self.zero,self.keep,self.abs,self.ins)
        CM_c=np.copy(self.mol.CM) # 先把原来的系数矩阵保留下来
        self.mol.datas['CM']=CM # 替换为投影后的系数矩阵
        # 计算各种性质
        if self.prop=='charge':
            values=self.atomCharge(self.chrg)
        if self.prop=='spin':
            values=self.atomSpin(self.chrg)
        # 替换为原本的系数矩阵
        self.mol.datas['CM']=CM_c 
        return values
    
    def atomCharge(self,chrg:str)->np.ndarray:
        """计算mulliken电子数"""
        caler=atomCharge.Calculator(self.mol)
        caler.chrg=chrg
        charges=caler.calculate()
        atomics=np.array(self.mol.atoms.atomics)
        eleNums=atomics-charges
        idxs=[a-1 for a in self.atoms]
        return eleNums[idxs]

    def atomSpin(self,chrg:str)->np.ndarray:
        """计算mulliken自旋"""
        caler=atomSpin.Calculator(self.mol)
        caler.chrg=chrg
        spins=caler.calculate()
        idxs=[a-1 for a in self.atoms]
        return spins[idxs]
    
    def atomEngs(self):
        """原子能量"""
        caler=atomEnergy.Calculator(self.mol)
        engs=caler.calculate()
        idxs=[a-1 for a in self.atoms]
        return engs[idxs]

    def freeValence(self):
        """计算原子自由价"""
        caler=freeValence.Calculator(self.mol)
        caler.calculate()

    def projVector(self)->np.ndarray:
        """轨道系数向量投影"""
        obts=self.mol.O_obts
        proj=np.zeros(shape=(len(self.atoms),len(obts),3)) #[原子/方向，轨道]
        for a,atm in enumerate(self.atoms):
            vect=self.vects[a]
            atom=self.mol.atom(atm)
            if atom.symbol=='H':
                proj[a,:]=0
                continue
            for o,obt in enumerate(obts):
                Cops=atom.get_pProj(vect,obt,False) # p轨道向量
                Cop=np.array(Cops).mean(axis=0)
                proj[a,o]=Cop
        return proj.sum(axis=1)
    
    def printRes(self):
        resStr=self.resStr()
        printer.info('方向电子自旋分布: ')
        printer.res(resStr)
        
    def resStr(self):
        values=self.calculate()
        return lutils.atomValueStr(self.mol,self.atoms,values)
