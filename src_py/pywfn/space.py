"""
空间工具

对应 rswfn 的 `space` 模块，提供格点生成与原子权重分配等空间相关计算。
"""

import numpy as np
from numpy.typing import NDArray

from pywfn import _core


def a2m_weits(
    iatm: int,
    atm_grids: NDArray[np.float64],
    atm_weits: NDArray[np.float64],
    atm_pos: NDArray[np.float64],
    atm_rad: NDArray[np.float64],
) -> NDArray[np.float64]:
    """将原子分区的格点权重重新分配到整个分子

    Args:
        iatm: 目标原子索引
        atm_grids: 原子分区的格点坐标
        atm_weits: 原子分区的格点权重
        atm_pos: 所有原子的位置
        atm_rad: 所有原子的半径

    Returns:
        NDArray[np.float64]: 重新分配后的权重
    """
    return _core.space.a2m_weits(iatm, atm_grids, atm_weits, atm_pos, atm_rad)  # type: ignore


def calc_EIEBA(
    fragA: NDArray[np.float64],
    fragB: NDArray[np.float64],
    rhoA: NDArray[np.float64],
    rhoB: NDArray[np.float64],
    gridsA: NDArray[np.float64],
    gridsB: NDArray[np.float64],
    weitsA: NDArray[np.float64],
    weitsB: NDArray[np.float64],
) -> tuple[float, float, float]:
    """计算两个片段之间的 EIEBA 值

    Args:
        fragA: 片段 A 的格点数据，每个元素为 [x, y, z, 值]
        fragB: 片段 B 的格点数据
        rhoA: 片段 A 的密度值
        rhoB: 片段 B 的密度值
        gridsA: 片段 A 的格点坐标
        gridsB: 片段 B 的格点坐标
        weitsA: 片段 A 的格点权重
        weitsB: 片段 B 的格点权重

    Returns:
        tuple[float, float, float]: 三个计算结果
    """
    return _core.space.calc_EIEBA(  # type: ignore
        fragA, fragB, rhoA, rhoB, gridsA, gridsB, weitsA, weitsB
    )


def get_grids(nx: int, ny: int, nz: int) -> NDArray[np.float64]:
    """生成三维格点坐标

    Args:
        nx: x 方向格点数
        ny: y 方向格点数
        nz: z 方向格点数

    Returns:
        NDArray[np.float64]: 形状为 (nx*ny*nz, 3) 的格点坐标
    """
    return np.array(_core.space.get_grids(nx, ny, nz))  # type: ignore
