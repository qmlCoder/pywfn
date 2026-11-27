"""
量子化学计算方法模块。

对应 rswfn 的 `qmcal` 模块，包含 Hartree-Fock（hf）、Kohn-Sham（ks）、
局部 π 电子 Kohn-Sham（lks）、半经验（se）、dHMO（dhmo/dhmo2）、
DIIS（diis）与初猜（guess）等子模块。
"""

from pywfn import _core
from pywfn.qmcal import dhmo, dhmo2, diis, guess, hf, ks, lks, se


def SlaterX(rho: float) -> float:
    """Slater 型原子轨道的指数换算函数

    按 Slater 规则把标度参数换算成对应的轨道指数，供 dHMO（双中心 Hückel）
    计算构造 Slater 型基函数时使用。

    Args:
        rho: 标度参数（无量纲）

    Returns:
        float: 换算后的 Slater 指数
    """
    return _core.qmcal.SlaterX(rho)  # type: ignore


__all__ = [
    "SlaterX",
    "dhmo",
    "dhmo2",
    "diis",
    "guess",
    "hf",
    "ks",
    "lks",
    "se",
]
