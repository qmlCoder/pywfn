"""
自洽场初猜

对应 rswfn 的 `qmcal::guess` 模块，提供原子轨道归属查询与叠加原子密度（SAD）
初猜密度矩阵的构造。
"""

import numpy as np
from numpy.typing import NDArray

from pywfn import _core
from pywfn.base.mole import Mole


def ao_owners(mole: Mole, form: str = "sph") -> list[int]:
    """返回每个原子轨道所属的原子索引

    Args:
        mole: 分子对象
        form: 基函数形式（"sph" 或 "car"）

    Returns:
        list[int]: 长度为原子轨道数的列表，每个元素是该轨道所属的原子索引
    """
    return _core.qmcal.guess.ao_owners(mole.core, form)


def use_sad_enabled() -> bool:
    """是否启用 SAD（叠加原子密度）初猜

    Returns:
        bool: 启用时返回 True，否则返回 False
    """
    return _core.qmcal.guess.use_sad_enabled()


def sad_dmat(
    smat: NDArray[np.float64],
    hmat: NDArray[np.float64],
    ao_atm: list[int],
    nucs: list[int],
    nele: int,
) -> NDArray[np.float64] | None:
    """用叠加原子密度（SAD）方法构造初始密度矩阵

    Args:
        smat: 重叠矩阵，形状为 (nato, nato)
        hmat: 核心哈密顿矩阵，形状为 (nato, nato)
        ao_atm: 每个原子轨道所属的原子索引
        nucs: 各原子的原子序数
        nele: 电子数

    Returns:
        NDArray[np.float64] | None: 初始密度矩阵；无法构造时返回 None
    """
    return _core.qmcal.guess.sad_dmat(smat, hmat, ao_atm, nucs, nele)  # type: ignore
