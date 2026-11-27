"""
数学工具

对应 rswfn 的 `maths` 模块，提供插值、矩阵运算、向量运算等基础功能。
"""

import numpy as np
from numpy.typing import NDArray

from pywfn import _core
from pywfn.base import Vec3


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


def arg_max(vals: list[float]) -> tuple[int, float]:
    """获取列表中的最大值及其索引

    Args:
        vals: 数值列表

    Returns:
        tuple[int, float]: (索引, 最大值)
    """
    return _core.maths.arg_max(vals)  # type: ignore


def arg_min(vals: list[float]) -> tuple[int, float]:
    """获取列表中的最小值及其索引

    Args:
        vals: 数值列表

    Returns:
        tuple[int, float]: (索引, 最小值)
    """
    return _core.maths.arg_min(vals)  # type: ignore


def bolzman_dist(ens: list[float], t: float) -> list[float]:
    """按玻尔兹曼分布计算各状态的概率

    Args:
        ens: 各状态的能量列表
        t: 温度

    Returns:
        list[float]: 各状态的概率
    """
    return _core.maths.bolzman_dist(ens, t)  # type: ignore


def calc_dist(p1: list[float], p2: list[float]) -> float:
    """计算两点之间的距离

    Args:
        p1: 第一个点坐标
        p2: 第二个点坐标
    """
    return _core.maths.calc_dist(p1, p2)  # type: ignore


def det(mat: NDArray[np.float64]) -> float:
    """计算矩阵的行列式"""
    return _core.maths.det(mat)  # type: ignore


def eigh(mat: NDArray[np.float64]) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """对对称矩阵做特征分解

    Args:
        mat: 对称矩阵

    Returns:
        tuple: (特征值数组, 特征向量矩阵)
    """
    vals, vecs = _core.maths.eigh(mat)  # type: ignore
    return vals, vecs


def get_plan_by_3points(p1: Vec3, p2: Vec3, p3: Vec3) -> list[float]:
    """由三个点确定一个平面

    Args:
        p1: 平面上的第一个点
        p2: 平面上的第二个点
        p3: 平面上的第三个点

    Returns:
        list[float]: 平面方程 ax + by + cz + d = 0 的系数 [a, b, c, d]
    """
    return _core.maths.get_plan_by_3points(p1.core, p2.core, p3.core)  # type: ignore


def inv(mat: NDArray[np.float64]) -> NDArray[np.float64]:
    """计算矩阵的逆矩阵"""
    return _core.maths.inv(mat)  # type: ignore


def sign(a: float, b: float) -> float:
    """返回与 b 同号的 a

    Args:
        a: 数值
        b: 决定符号的数值
    """
    return _core.maths.sign(a, b)  # type: ignore


def vec_cross(vec1: list[float], vec2: list[float]) -> list[float]:
    """计算两个向量的叉积"""
    return _core.maths.vec_cross(vec1, vec2)  # type: ignore


def vec_dist(vec1: list[float], vec2: list[float]) -> float:
    """计算两个向量的距离"""
    return _core.maths.vec_dist(vec1, vec2)  # type: ignore


def vec_len(vec: list[float]) -> float:
    """计算向量的长度"""
    return _core.maths.vec_len(vec)  # type: ignore


def vec_norm(vec: list[float]) -> list[float]:
    """将向量归一化"""
    return _core.maths.vec_norm(vec)  # type: ignore


def vec_project(vec1: list[float], vec2: list[float]) -> list[float]:
    """计算向量在另一个向量方向上的投影"""
    return _core.maths.vec_project(vec1, vec2)  # type: ignore


def vec_psd(nums: list[float]) -> float:
    """计算向量的平方和"""
    return _core.maths.vec_psd(nums)  # type: ignore
