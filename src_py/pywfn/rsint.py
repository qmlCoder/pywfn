import numpy as np

from pywfn import _core
from pywfn.base import Mole, Stm


class RsInt:
    def __init__(self, mole: Mole) -> None:
        self.mole = mole
        self.core = _core.rsint.RsInt(mole.core)  # type: ignore

    def smat(self, form: str) -> np.ndarray:
        """计算重叠矩阵"""
        return self.core.smat(form)

    def tmat(self, form: str) -> np.ndarray:
        """计算动能矩阵"""
        return self.core.tmat(form)

    def vmat(self, form: str) -> np.ndarray:
        """计算势能矩阵"""
        return self.core.vmat(form)

    def mat1e(self, cals: list[int], form: str) -> np.ndarray:
        """计算单电子积分矩阵

        Args:
            cals: 指定需要计算的矩阵，[calc_s, calc_t, calc_v]，1 表示计算，0 表示跳过
            form: 积分格式，"sph" 或 "car"

        Returns:
            三个二维数组 [S, T, V]
        """
        return self.core.mat1e(cals, form)

    def mat2e(self, form: str) -> np.ndarray:
        """计算双电子积分矩阵

        Args:
            form: 积分格式，"sph" 或 "car"
        """
        return self.core.mat2e(form)

    def mat1e_lcs(self, stms: dict[int, Stm], cals: list[int]) -> np.ndarray:
        """计算局部坐标系下的单电子积分矩阵"""
        stms = {atm: stm.core for atm, stm in stms.items()}
        return self.core.mat1e_lcs(stms, cals)
