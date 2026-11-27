"""
计算分子轨道空间波函数

基函数的波函数、原子轨道的波函数、原子的波函数

镜像自 Rust 模块 `src_rs/gridprop/wfnfunc.rs` 的结构，不包含实际代码逻辑。
"""

import numpy as np

from pywfn import _core
from pywfn.base import Mole

Grids = list[list[float]] | np.ndarray


class Calculator:
    def __init__(self, mole: Mole) -> None:
        self.mole = mole
        self.atms: list[int] = []
        self.core = _core.props.gridprop.wfnfunc.Calculator(mole.core)  # type: ignore

    def ato_wfn(
        self, grids: Grids, level: int
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """计算所有原子轨道在给定格点上的波函数值。

        Args:
            grids: 三维坐标点列表 [[x, y, z], ...]
            level: 计算级别（0=仅波函数值，1=值+梯度，2=值+梯度+Hessian）

        Returns:
            每个格点对应一个元素，包含该点处每个原子轨道的
            (波函数值, 梯度, Hessian) 向量
        """
        return self.core.ato_wfn(grids, level)

    def obt_wfn(
        self, obt: int, grids: Grids, level: int
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """计算指定分子轨道在给定格点上的波函数值。

        Args:
            obt: 分子轨道索引（从 0 开始）
            grids: 三维坐标点列表 [[x, y, z], ...]
            level: 计算级别（0=仅波函数值，1=值+梯度，2=值+梯度+Hessian）

        Returns:
            (波函数值, 梯度, Hessian) 三元组，每个元素为对应格点的列表
        """
        return self.core.obt_wfn(obt, grids, level)


def gtf(
    xyzs: Grids,
    lmn: list[int],
    alp: float,
    level: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """计算高斯型轨道（GTO）在给定格点上的值。

    Args:
        xyzs: 以原子为中心的笛卡尔坐标点列表
        lmn: 角动量指数 [l, m, n]，对应 x^l * y^m * z^n
        alp: 高斯指数参数
        level: 计算级别（0=仅波函数值，1=值+梯度，2=值+梯度+Hessian）

    Returns:
        (波函数值, 梯度, Hessian) 三元组
    """
    return _core.props.gridprop.wfnfunc.gtf(xyzs, lmn, alp, level)  # type: ignore


def cgf(
    xyzs: Grids,
    lmn: list[int],
    alps: list[float],
    coes: list[float],
    level: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """计算收缩高斯函数（CGF）在给定格点上的值。

    CGF 是多个 GTO 的线性组合，用于表示原子轨道。

    Args:
        xyzs: 以原子为中心的笛卡尔坐标点列表
        lmn: 角动量指数 [l, m, n]
        alps: 高斯指数参数列表
        coes: 收缩系数列表
        level: 计算级别（0=仅波函数值，1=值+梯度，2=值+梯度+Hessian）

    Returns:
        每个格点上的 (波函数值, 梯度, Hessian) 向量列表
    """
    return _core.props.gridprop.wfnfunc.cgf(xyzs, lmn, alps, coes, level)  # type: ignore

def ato_wfn(mole: Mole, grid: list[float], level: int = 0):
    """计算单个原子在格点上的波函数值

    Args:
        mole: 分子对象
        grid: 格点坐标
        level: 计算级别（0=仅值，1=值+梯度，2=值+梯度+Hessian）

    Returns:
        (函数值数组, 梯度数组, Hessian 数组) 三元组
    """
    return _core.props.gridprop.wfnfunc.ato_wfn(mole.core, grid, level)  # type: ignore


def stf(xyz: list[float], nuc: int, Z: float, nlm: tuple[int, int, int]) -> float:
    """计算 Slater 型轨道的值

    Args:
        xyz: 格点坐标
        nuc: 原子序数
        Z: 有效核电荷
        nlm: 量子数 (n, l, m)

    Returns:
        float: Slater 型轨道在格点上的值
    """
    return _core.props.gridprop.wfnfunc.stf(xyz, nuc, Z, nlm)  # type: ignore