"""
Rys 求积的根与权重

对应 rswfn 的 `datas::roots` 模块。
"""

import numpy as np

from pywfn import _core

# Rys 多项式求积的 alpha 参数
RYS_ALPS: list[float] = _core.datas.roots.RYS_ALPS()  # type: ignore

# 各阶 Rys 求积的根
RYS_ROOTS_1: list[list[float]] = _core.datas.roots.RYS_ROOTS_1()  # type: ignore
RYS_ROOTS_2: list[list[float]] = _core.datas.roots.RYS_ROOTS_2()  # type: ignore
RYS_ROOTS_3: list[list[float]] = _core.datas.roots.RYS_ROOTS_3()  # type: ignore
RYS_ROOTS_4: list[list[float]] = _core.datas.roots.RYS_ROOTS_4()  # type: ignore
RYS_ROOTS_5: list[list[float]] = _core.datas.roots.RYS_ROOTS_5()  # type: ignore
RYS_ROOTS_6: list[list[float]] = _core.datas.roots.RYS_ROOTS_6()  # type: ignore
RYS_ROOTS_7: list[list[float]] = _core.datas.roots.RYS_ROOTS_7()  # type: ignore
RYS_ROOTS_8: list[list[float]] = _core.datas.roots.RYS_ROOTS_8()  # type: ignore
RYS_ROOTS_9: list[list[float]] = _core.datas.roots.RYS_ROOTS_9()  # type: ignore

# 各阶 Rys 求积的权重
RYS_WEITS_1: list[list[float]] = _core.datas.roots.RYS_WEITS_1()  # type: ignore
RYS_WEITS_2: list[list[float]] = _core.datas.roots.RYS_WEITS_2()  # type: ignore
RYS_WEITS_3: list[list[float]] = _core.datas.roots.RYS_WEITS_3()  # type: ignore
RYS_WEITS_4: list[list[float]] = _core.datas.roots.RYS_WEITS_4()  # type: ignore
RYS_WEITS_5: list[list[float]] = _core.datas.roots.RYS_WEITS_5()  # type: ignore
RYS_WEITS_6: list[list[float]] = _core.datas.roots.RYS_WEITS_6()  # type: ignore
RYS_WEITS_7: list[list[float]] = _core.datas.roots.RYS_WEITS_7()  # type: ignore
RYS_WEITS_8: list[list[float]] = _core.datas.roots.RYS_WEITS_8()  # type: ignore
RYS_WEITS_9: list[list[float]] = _core.datas.roots.RYS_WEITS_9()  # type: ignore


def dub_roots() -> None:
    """打印 Rys 求积数据表"""
    _core.datas.roots.dub_roots()  # type: ignore


def hermite(n: int) -> np.ndarray:
    """计算 Hermite 求积的根和权重

    Args:
        n: 求积点数

    Returns:
        np.ndarray: 形状为 (n, 2) 的数组，每行为 (根, 权重)
    """
    return np.array(_core.datas.roots.hermite(n))  # type: ignore


def legendre(n: int) -> np.ndarray:
    """计算 Legendre 求积的根和权重

    Args:
        n: 求积点数

    Returns:
        np.ndarray: 形状为 (n, 2) 的数组，每行为 (根, 权重)
    """
    return np.array(_core.datas.roots.legendre(n))  # type: ignore


def rys_roots(n: int, alp: float) -> np.ndarray:
    """计算 Rys 求积的根和权重

    Args:
        n: 求积点数
        alp: Rys 参数

    Returns:
        np.ndarray: 形状为 (n, 2) 的数组，每行为 (根, 权重)
    """
    return np.array(_core.datas.roots.rys_roots(n, alp))  # type: ignore
