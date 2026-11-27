from pywfn import _core
from pywfn.base.mole import Mole
from pywfn.base import Stm


def get_integ_x(zi: int, zj: int, pi: int, pj: int, dist: float) -> float:
    """计算两个位于 x 轴上的基函数的重叠积分（距离 > 3.0 时返回 0）"""
    return _core.qmcal.dhmo2.get_integ_x(zi, zj, pi, pj, dist)


def get_integ_d(
    zi: int, zj: int, di: "list[float]", dj: "list[float]", stm: Stm, dist: float
) -> float:
    """计算两个任意朝向的 Pz 轨道的重叠积分"""
    return _core.qmcal.dhmo2.get_integ_d(zi, zj, di, dj, stm.core, dist)


def dhmo_calc_v2(mole: Mole) -> None:
    """dHMO 方法版本 2，计算分子中的所有原子轨道（输出到日志）"""
    _core.qmcal.dhmo2.dhmo_calc_v2(mole.core)
