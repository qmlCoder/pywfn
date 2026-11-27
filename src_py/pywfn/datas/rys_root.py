"""
Rys 多项式求根的批量数据

对应 rswfn 的 `datas::rys_root` 模块。
"""

import numpy as np

from pywfn import _core

# 大 X 区间的根与权重
POLY_LARGEX_RT: list[float] = _core.datas.rys_root.POLY_LARGEX_RT()  # type: ignore
POLY_LARGEX_WW: list[float] = _core.datas.rys_root.POLY_LARGEX_WW()  # type: ignore

# 小 X 区间的根与权重
POLY_SMALLX_R0: list[float] = _core.datas.rys_root.POLY_SMALLX_R0()  # type: ignore
POLY_SMALLX_R1: list[float] = _core.datas.rys_root.POLY_SMALLX_R1()  # type: ignore
POLY_SMALLX_W0: list[float] = _core.datas.rys_root.POLY_SMALLX_W0()  # type: ignore
POLY_SMALLX_W1: list[float] = _core.datas.rys_root.POLY_SMALLX_W1()  # type: ignore


def rys_root1(x: float) -> tuple[list[float], list[float]]:
    """计算 Rys 多项式的一阶根和权重"""
    return _core.datas.rys_root.rys_root1(x)  # type: ignore


def rys_root2(x: float) -> tuple[list[float], list[float]]:
    """计算 Rys 多项式的二阶根和权重"""
    return _core.datas.rys_root.rys_root2(x)  # type: ignore


def rys_root3(x: float) -> tuple[list[float], list[float]]:
    """计算 Rys 多项式的三阶根和权重"""
    return _core.datas.rys_root.rys_root3(x)  # type: ignore


def rys_root4(x: float) -> tuple[list[float], list[float]]:
    """计算 Rys 多项式的四阶根和权重"""
    return _core.datas.rys_root.rys_root4(x)  # type: ignore


def rys_root5(x: float) -> tuple[list[float], list[float]]:
    """计算 Rys 多项式的五阶根和权重"""
    return _core.datas.rys_root.rys_root5(x)  # type: ignore


def rys_roots(n: int, x: float) -> np.ndarray:
    """计算任意阶 Rys 多项式的根和权重

    Args:
        n: 阶数
        x: Rys 多项式参数

    Returns:
        np.ndarray: 形状为 (n, 2) 的数组，每行为 (根, 权重)
    """
    return np.array(_core.datas.rys_root.rys_roots(n, x))  # type: ignore
