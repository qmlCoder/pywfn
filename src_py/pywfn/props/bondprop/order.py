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
        self.core = _core.props.bondprop.order.Calculator(mole.core)  # type: ignore # 核心计算器

    def mayer(self) -> np.ndarray:
        """计算Mayer键级矩阵"""
        return self.core.mayer()

    def wiberg(self) -> np.ndarray:
        """计算Wiberg键级矩阵"""
        return self.core.wiberg()

    def order_obt(self, otype: str, iobt: int) -> np.ndarray:
        """计算OBT键级矩阵"""
        return self.core.order_obt(otype, iobt)

    def pocv(
        self,
        dirs: dict[int, list[float]],
        keep_other_atm: bool,
        keep_other_sym: bool,
        otype: str = "mayer",
    ) -> np.ndarray:
        """计算带有投影方向的键级

        Args:
            dirs: 方向字典, {原子索引: 方向向量}
            keep_other_atm: 是否保留其他原子
            keep_other_sym: 是否保留其他对称性
            otype: 键级类型, 默认为 "mayer"

        Returns:
            np.ndarray: 原子间键级矩阵, 形状为 [natm, natm]
        """
        return self.core.pocv(dirs, keep_other_atm, keep_other_sym, otype)  # type: ignore

    def mocv(
        self, stms: dict[int, Stm], atos: list[int], otype: str = "mayer"
    ) -> np.ndarray:
        stms = {atm: stm.core for atm, stm in stms.items()}
        return self.core.mocv(stms, atos, otype)

    def bound(
        self,
        atm: int,
        dir: list[float],
        nebs: list[int],
        otype: str = "mayer",
        mtype: str = "pocv",
    ) -> list[float]:
        """计算与指定原子相邻的键的束缚键级

        Args:
            atm: 指定原子编号
            dir: 原子投影的方向 [3] (x,y,z)
            nebs: 相邻原子列表
            otype: 键级类型, 默认为 "mayer"
            mtype: 分解方法, 可选 "pocv" 或 "mocv", 默认为 "pocv"

        Returns:
            list[float]: 与 `nebs` 一一对应的束缚键级
        """
        return self.core.bound(atm, dir, nebs, otype, mtype)

    def pi_pocv(
        self, otype: str = "mayer"
    ) -> tuple[dict[int, list[float]], np.ndarray]:
        """使用pocv方法计算pi键级矩阵"""
        dirs, omat = self.core.pi_pocv(otype)
        return dirs, omat

    def pi_mocv(
        self, otype: str = "mayer"
    ) -> tuple[dict[int, Stm], list[int], np.ndarray]:
        """使用mocv方法计算pi键级"""
        stms, atos, omat = self.core.pi_mocv(otype)
        stms = {atm: Stm(stm) for atm, stm in stms.items()}
        return stms, atos, omat  # type: ignore

    def decom(self, ai: int, aj: int, otype: str = "mayer") -> list[float]:
        return self.core.decom(ai, aj, otype)


def mayer(dmat: np.ndarray, smat: np.ndarray) -> np.ndarray:
    """计算Mayer键级

    Args:
        dmat: 密度矩阵
        smat: 重叠矩阵

    Returns:
        np.ndarray: Mayer键级矩阵
    """
    return _core.props.bondprop.order.mayer(dmat, smat)  # type: ignore


def wiberg(dmat: np.ndarray, smat: np.ndarray) -> np.ndarray:
    """计算Wiberg键级

    Args:
        dmat: 密度矩阵
        smat: 重叠矩阵

    Returns:
        np.ndarray: Wiberg键级矩阵
    """
    return _core.props.bondprop.order.wiberg(dmat, smat)  # type: ignore


def ato2atm(omat: np.ndarray, uplo: list[list[int]]) -> np.ndarray:
    """将原子轨道键级矩阵转换为原子间键级矩阵

    Args:
        omat: 原子轨道键级矩阵
        uplo: 原子轨道索引

    Returns:
        np.ndarray: 原子间键级矩阵
    """
    return _core.props.bondprop.order.ato2atm(omat, uplo)  # type: ignore
