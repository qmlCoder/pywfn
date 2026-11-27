"""
分子轨道系数变换

对应 rswfn 的 `props::orbtprop::mocv` 模块。
"""

import numpy as np

from pywfn import _core
from pywfn.base import Mole

core = _core.props.orbtprop.mocv  # type: ignore


class OrderDType:
    """键级分解的类型"""

    SI = _core.props.orbtprop.mocv.OrderDType.SI  # type: ignore # 总键级
    PI = _core.props.orbtprop.mocv.OrderDType.PI  # type: ignore # π 键级
    DE = _core.props.orbtprop.mocv.OrderDType.DE  # type: ignore # σ 键级
    PH = _core.props.orbtprop.mocv.OrderDType.PH  # type: ignore # 相位键级


class MocvDest:
    """分子轨道系数的用途"""

    @staticmethod
    def PiProperty():
        """用于计算 π 电子相关的性质

        Returns:
            MocvDest: 用途对象
        """
        return _core.props.orbtprop.mocv.MocvDest.PiProperty()  # type: ignore

    @staticmethod
    def OrderDecom(dtype: OrderDType):
        """用于分解键级

        Args:
            dtype: 需要分解的键级类型

        Returns:
            MocvDest: 用途对象
        """
        return _core.props.orbtprop.mocv.MocvDest.OrderDecom(dtype)  # type: ignore

    @staticmethod
    def RoundOrder():
        """用于取整键级

        Returns:
            MocvDest: 用途对象
        """
        return _core.props.orbtprop.mocv.MocvDest.RoundOrder()  # type: ignore


def mocv(L: np.ndarray, rcoefs: np.ndarray, keeps: np.ndarray) -> np.ndarray:
    """计算分子轨道系数变换矩阵

    Args:
        L: 分子轨道系数矩阵
        rcoefs: 参考系数矩阵
        keeps: 保留的轨道掩码

    Returns:
        np.ndarray: 变换后的系数矩阵
    """
    return core.mocv(L, rcoefs, keeps)  # type: ignore


def mocv_atos(mole: Mole, dest) -> list[int]:
    """获取参与分子轨道系数变换的原子索引

    Args:
        mole: 分子对象
        dest: 分子轨道系数的用途

    Returns:
        list[int]: 原子索引列表
    """
    return core.mocv_atos(mole.core, dest)  # type: ignore


def mocv_syms(dest, geom: int) -> list[str]:
    """获取参与分子轨道系数变换的基函数符号

    Args:
        dest: 分子轨道系数的用途
        geom: 几何结构类型

    Returns:
        list[str]: 基函数符号列表
    """
    return core.mocv_syms(dest, geom)  # type: ignore


def get_tmat(num: int, L: np.ndarray) -> np.ndarray:
    """获取变换矩阵

    Args:
        num: 原子序数
        L: 局部坐标系

    Returns:
        np.ndarray: 变换矩阵
    """
    return core.get_tmat(num, L)  # type: ignore


def get_tmat_sph(num: int, L: np.ndarray) -> np.ndarray:
    """获取球谐变换矩阵"""
    return core.get_tmat_sph(num, L)  # type: ignore


def get_tmat_car(num: int, L: np.ndarray) -> np.ndarray:
    """获取笛卡尔变换矩阵"""
    return core.get_tmat_car(num, L)  # type: ignore


def tmat_s(L: np.ndarray) -> np.ndarray:
    """获取s轨道的变换矩阵"""
    return core.tmat_s(L)  # type: ignore


def tmat_p(L: np.ndarray) -> np.ndarray:
    """获取p轨道的变换矩阵"""
    return core.tmat_p(L)  # type: ignore


def tmat_d(L: np.ndarray) -> np.ndarray:
    """获取d轨道的变换矩阵"""
    return core.tmat_d(L)  # type: ignore


def tmat_f(L: np.ndarray) -> np.ndarray:
    """获取f轨道的变换矩阵"""
    return core.tmat_f(L)  # type: ignore


def tmat_g(L: np.ndarray) -> np.ndarray:
    """获取g轨道的变换矩阵"""
    return core.tmat_g(L)  # type: ignore


def get_nmat(lmns: list[list[int]]) -> np.ndarray:
    """获取n矩阵"""
    return core.get_nmat(lmns)  # type: ignore


def nmat_d() -> np.ndarray:
    """获取d轨道的n矩阵"""
    return core.nmat_d()  # type: ignore


def nmat_f() -> np.ndarray:
    """获取f轨道的n矩阵"""
    return core.nmat_f()  # type: ignore


def nmat_g() -> np.ndarray:
    """获取g轨道的n矩阵"""
    return core.nmat_g()  # type: ignore
