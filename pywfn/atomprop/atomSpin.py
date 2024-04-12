"""
此脚本用来计算mulliken自旋
"""
import numpy as np
from pywfn.base import Mol,Atom

from pywfn.utils import printer
from pywfn.atomprop import atomCharge, lutils,AtomCaler
from typing import Literal


class Calculator(AtomCaler):
    def __init__(self,mol:"Mol") -> None:
        self.mol=mol
        self.logTip='mulliken 电子自旋分布:'
        self.chrg:Literal['mulliken','lowdin']='mulliken'
    
    def get_Es(self,obts:list[int])->list[float]:
        elects=lutils.get_ects(self.mol,obts,self.mol.CM)
        return elects

    def calculate(self)->np.ndarray:
        """计算所有原子的自旋"""
        assert not self.mol.isOpenShell,'非开壳层分子无法计算自旋'
        obtNum=self.mol.CM.shape[0] # 系数矩阵行数，基函数数量
        CM_a=self.mol.CM[:,:obtNum]
        CM_b=self.mol.CM[:,obtNum:]
        obtEcts=self.mol.obtEcts
        a_obt=[i for i,e in enumerate(obtEcts) if (e!=0 and i< obtNum)]
        b_obt=[i for i,e in enumerate(obtEcts) if (e!=0 and i>=obtNum)]
        # 记录原本的分子属性
        CM_r=self.mol.CM.copy()
        obt_r=self.mol.O_obts
        # 通过临时修改分子属性，计算想要结果
        if self.chrg=='mulliken':
            caler=atomCharge.Calculator(self.mol)
            caler.chrg='mulliken'
        if self.chrg=='lowdin':
            caler=atomCharge.Calculator(self.mol)
            caler.chrg='lowdin'
        # 将长方形的系数矩阵分为两个正方形分别计算
        self.mol.datas['CM']=CM_a.copy()
        self.mol.datas['O_obts']=a_obt
        a_Ects=caler.calculate()
        self.mol.datas['CM']=CM_b.copy()
        self.mol.datas['O_obts']=b_obt
        b_Ects=caler.calculate()
        # 恢复分子属性
        self.mol.datas['CM']=CM_r
        self.mol.datas['O_obts']=obt_r

        # a_Ects=np.array(self.get_Es(a_obt))
        # b_Ects=np.array(self.get_Es(b_obt))
        return a_Ects-b_Ects
        
    def resStr(self):
        elects=self.calculate()
        atoms=lutils.atomIdxs(self.mol.atoms)
        return lutils.atomValueStr(self.mol,atoms,elects)