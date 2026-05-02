"""
计算分子空间电子密度
整个分子只有网格点，但是没有权重
对于每个原子，根据自身的坐标及权重插值到整个空间的坐标及权重
？ 波函数值是否也要插值呢？

空间中一个格点的权重是否应该对所有分子都一致？
"""
from pywfn.base.mole import Mole
import numpy as np
from pywfn import _core


class Calculator:
    def __init__(self, mole: Mole) -> None:
        self.caler = _core.gridprop.density.Calculator(mole.core)  # type: ignore # 核心计算器

    def mol_rho_cm(self, grids: np.ndarray, level: int):  # 以分子轨道的方式计算电子密度，可以计算每个分子轨道的电子密度
        """计算分子电子密度"""
        return self.caler.mol_rho_cm(grids, level)
    
    def mol_rho_dm(self,grids:np.ndarray,level:int):
        return self.caler.mol_rho_dm(grids,level)

    def pi_rho(self, grids: np.ndarray):
        return self.caler.pi_rho(grids)

    def ato_rho(self, grids: np.ndarray, level: int):
        """计算所有原子轨道的电子密度"""
        return self.caler.ato_rho(grids, level)
