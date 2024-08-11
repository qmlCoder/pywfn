"""
此脚本用来计算mulliken自旋
自旋是α电子-β电子
可以通过改变轨道的占据情况来分别计算α和β电子数
"""

from pywfn.base import Mol,Atom
from pywfn.utils import printer
from pywfn.atomprop import charge, lutils,AtomCaler
from pywfn.atomprop.charge import Chrgs
from pywfn.maths import CM2PM

from typing import Literal
import numpy as np

class Calculator(AtomCaler):
    def __init__(self,mol:"Mol") -> None:
        self.mol=mol
        self.logTip='mulliken 电子自旋分布:'

    def calculate(self,chrg:Chrgs='mulliken')->np.ndarray:
        """计算所有原子的自旋"""
        assert self.mol.open,'非开壳层分子无法计算自旋'
        obtNum=self.mol.CM.shape[0] # 系数矩阵行数，基函数数量
        occs_old=self.mol.obtOccs # 记录原本的占据情况
        a_occs=occs_old.copy() # 当你需要修改一个变量的时候，
        b_occs=occs_old.copy()
        a_occs[obtNum:]=[False]*obtNum
        b_occs[:obtNum]=[False]*obtNum

        caler=charge.Calculator(self.mol)

        # 将长方形的系数矩阵分为两个正方形分别计算
        PMa=CM2PM(self.mol.CM.copy(),a_occs,1)
        PMb=CM2PM(self.mol.CM.copy(),b_occs,1)
        a_Ects=caler.charge(chrg=chrg,PM=PMa)
        b_Ects=caler.charge(chrg=chrg,PM=PMb)
        # 恢复分子属性
        return -(a_Ects-b_Ects)
    
    def onShell(self):
        while True:
            printer.options('原子自旋',{
                '1':'mulliken电子自旋',
                '2':'lowdin电子自旋',
            })
            opt=input('请输入自旋类型:')
            if opt=='1':
                spins=self.calculate(chrg='mulliken')
                for i,spin in enumerate(spins):
                    printer.res(f'{i+1}: {spin}')
            elif opt=='2':
                spins=self.calculate(chrg='lowdin')
                for i,spin in enumerate(spins):
                    printer.res(f'{i+1}: {spin}')
            else:
                break