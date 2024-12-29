"""
计算轨道布局
"""
import numpy as np

from pywfn.base import Mol
from pywfn.maths import flib
from pywfn.orbtprop import decom

class Calculator:
    def __init__(self,mol:Mol) -> None:
        self.mol=mol

    def eleMat(self)->np.ndarray: # 电子布局矩阵
        """计算与分子轨道系数矩阵对应的电子分布矩阵"""
        # 使用法向量可以计算每个分子的pi电子分布
        
        obts=self.mol.O_obts

        nobt=len(obts)
        CM=self.mol.CM.copy()
        nmat,nobt=CM.shape
        NM=flib.eleMat(nmat,nobt,CM,self.mol.SM)*self.mol.oE
        return NM
    
    def piEleMatDecom(self): # 分解出pi分子轨道
        """计算与分子轨道系数矩阵对应的电子分布矩阵"""
        CMo=self.mol.CM.copy()
        CMt=decom.Calculator(self.mol).pi_decom()
        self.mol.props.set('CM',CMt)
        NM=self.eleMat()
        self.mol.props.set('CM',CMo)
        return NM