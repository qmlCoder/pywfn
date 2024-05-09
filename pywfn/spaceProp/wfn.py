"""
计算空间波函数
"""
from pywfn.base import Mol
from pywfn.maths import Gto
import numpy as np

class Calculator:
    def __init__(self) -> None:
        self.mol:Mol=None

    def molWfn(self,pos:np.ndarray):
        """计算分子波函数"""
        pass