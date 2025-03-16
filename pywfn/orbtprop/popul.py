"""
计算轨道布局
"""
import numpy as np

from pywfn.base import Mol

from pywfn.orbtprop import decom

class Calculator:
    def __init__(self,mol:Mol) -> None:
        self.mol=mol

    def eleMat(self,virtual:bool=True)->np.ndarray: # 电子布局矩阵
        """计算电子布局矩阵

        Args:
            virtual (bool, optional): 是否包含空轨道. Defaults to False.

        Returns:
            np.ndarray: 每个基函数在每个分子轨道内的占据数
        """
        # 使用法向量可以计算每个分子的pi电子分布
        from pywfn.maths import rlib
        if virtual:
            CM=self.mol.CM.copy()
        else:
            obts=self.mol.O_obts
            CM=self.mol.CM[:,obts].copy()
        
        # nmat,nobt=CM.shape
        # NM=flib.eleMat(nmat,nobt,CM,self.mol.SM)*self.mol.oE
        # return NM
        NM=rlib.ele_mat_rs(CM,self.mol.SM) # type: ignore
        return np.array(NM)*self.mol.oE
    
    def piEleMatDecom(self,dtype:str,virtual:bool=True)->np.ndarray: # 分解出pi分子轨道
        """计算轨道分解法计算出来的pi电子布局矩阵

        Args:
            dtype (str): 轨道分解计算pi轨道的类型，atom/bond
            virtual (bool, optional): 是否包含空轨道. Defaults to False.

        Returns:
            np.ndarray: 每个基函数的pi电子成分在每个分子轨道内的占据数
        """
        CMo=self.mol.CM.copy()
        CMt=decom.Calculator(self.mol).pi_decom(dtype)
        self.mol.coefs._CM=CMt
        NM=self.eleMat(virtual=virtual)
        self.mol.coefs._CM=CMo
        return NM