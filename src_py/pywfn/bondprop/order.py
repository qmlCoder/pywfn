"""
键级也不止一种，都在这里实现吧
"""

import numpy as np

from pywfn import _core
from pywfn.base import Stm
from pywfn.base.mole import Mole


class Calculator:
    def __init__(self, mole: Mole) -> None:
        self.mole = mole
        self.core = _core.bondprop.order.Calculator(mole.core)  # type: ignore # 核心计算器

    # mayer键级
    def mayer(self):
        """计算mayer键级"""
        omat = self.core.mayer()
        return omat

    def wiberg(self):
        omat = self.core.wiberg()
        return omat

    def pocv(
        self, dirs: dict[int, list[float]], keep_other_atm: bool, keep_other_sym: bool
    ) -> float:
        """计算带有方向的Mayer键级，指定一个键的多个方向

        Args:
            bonds (list[int]): 指定要计算哪个键级

        Returns:
            np.ndarray: 返回数组形状为:[d,6](a1,a2,x,y,z,v),其中d为键级方向数量
        """
        return self.caler.pocv(dirs, keep_other_atm, keep_other_sym, "mayer")  # type: ignore

    def mocv(self, stms: dict[int, Stm], ctype: str = "mulliken"):
        """分解键级"""
        stms = {atm: stm.core for atm, stm in stms.items()}
        return self.core.mocv(stms, ctype)

    def bound(
        self, atm: int, dir: list[float], nebs: list[int], otype: str = "mayer"
    ) -> tuple[list[tuple[int, int]], np.ndarray]:
        """计算与指定原子相邻的键的束缚键级
        Args:
            atm (int): 指定原子编号
            dir (list[float]): 指定原子投影的方向，形状为[3,](x,y,z)
        """
        return self.core.bound(atm, dir, nebs, otype)

    def pi_pocv(
        self, otype: str = "mayer"
    ) -> tuple[dict[int, list[float]], np.ndarray]:
        """使用pocv方法计算pi键级矩阵"""
        dirs, omat = self.core.pi_pocv(otype)
        return dirs, omat

    def pi_mocv(self, otype: str = "mayer"):
        """使用mocv方法计计算pi键级"""
        stms, atos, omat = self.core.pi_mocv(otype)
        stms = {atm: Stm(stm) for atm, stm in stms.items()}
        return stms, atos, omat  # type: ignore


def mayer(dmat, smat):
    return _core.bondprop.order.mayer(dmat, smat)  # type: ignore


def wiberg(dmat, smat):
    return _core.bondprop.order.wiberg(dmat, smat)  # type: ignore


def ato2atm(omat, uplo):
    return _core.bondprop.order.ato2atm(omat, uplo)  # type: ignore
