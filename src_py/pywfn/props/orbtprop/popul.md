"""
计算轨道布局
"""
import numpy as np

from pywfn.base.mole import Mole

from pywfn.orbtprop import decom

class Calculator:
    def __init__(self,mole:Mole) -> None:
        self.mole=mole

    def eleMat(self,virtual:bool=True)->np.ndarray: # 电子布局矩阵
        """计算电子布局矩阵

        Args:
            virtual (bool, optional): 是否包含空轨道. Defaults to False.

        Returns:
            np.ndarray: 每个基函数在每个分子轨道内的占据数
        """
        # 使用法向量可以计算每个分子的pi电子分布
        from pywfn import core
        if virtual:
            CM=self.mole.CM.copy()
        else:
            obts=self.mole.O_obts
            CM=self.mole.CM[:,obts].copy()
        NM=core.matrix.ele_mat(CM,self.mol.SM) # type: ignore
        return np.array(NM)*self.mol.oE # type: ignore
    
    def piEleMatDecom(self,dtype:str,virtual:bool=True)->np.ndarray: # 分解出pi分子轨道
        """计算轨道分解法计算出来的pi电子布局矩阵

        Args:
            dtype (str): 轨道分解计算pi轨道的类型，atom/bond
            virtual (bool, optional): 是否包含空轨道. Defaults to False.

        Returns:
            np.ndarray: 每个基函数的pi电子成分在每个分子轨道内的占据数
        """
        CMo=self.mole.CM.copy()
        CMt=decom.Calculator(self.mole).pi_decom(dtype)
        self.mole.coefs._CM=CMt
        NM=self.eleMat(virtual=virtual)
        self.mole.coefs._CM=CMo
        return NM