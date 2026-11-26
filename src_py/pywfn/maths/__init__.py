import numpy as np
from numpy.typing import NDArray

from pywfn import _core


def lag_intpol(xs: list[float], ys: list[float], ts: list[float]) -> list[float]:
    """拉格朗日插值

    Args:
        xs: 已知点的 x 坐标列表
        ys: 已知点的 y 坐标列表
        ts: 需要插值的目标点列表

    Returns:
        插值结果列表
    """
    return _core.maths.lag_intpol(xs, ys, ts)  # type: ignore


def calc_dmat(
    cmat: NDArray[np.float64], obt_occs: NDArray[np.float64]
) -> NDArray[np.float64]:
    """计算密度矩阵

    Args:
        cmat: 系数矩阵
        obt_occs: 轨道占据数

    Returns:
        密度矩阵
    """
    return _core.maths.calc_dmat(cmat, obt_occs)  # type: ignore


def matrix_nhalf(mat: NDArray[np.float64]) -> NDArray[np.float64]:
    """计算矩阵的负一半次幂"""
    return _core.maths.matrix_nhalf(mat)  # type: ignore


def matrix_phalf(mat: NDArray[np.float64]) -> NDArray[np.float64]:
    """计算矩阵的正一半次幂"""
    return _core.maths.matrix_phalf(mat)  # type: ignore


def fac(n: int) -> int:
    """计算阶乘 n!"""
    return _core.maths.fac(n)  # type: ignore


def fac2(n: int) -> int:
    """计算双阶乘 n!!"""
    return _core.maths.fac2(n)  # type: ignore
