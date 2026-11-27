from dataclasses import dataclass

from pywfn import _core
from pywfn.base.mole import Mole


@dataclass
class PiEnergy:
    """π/σ/πσ 能量分解结果（闭壳层，杂化泛函）

    满足 ``e_pi + e_sigma + e_cross == e_elec``。
    """

    e_pi: float
    """π 自能量"""
    e_sigma: float
    """σ 自能量"""
    e_cross: float
    """π-σ 交叉项"""
    e_cross_coul: float
    """交叉项的库仑部分"""
    e_cross_exch: float
    """交叉项的精确交换部分"""
    e_cross_xc: float
    """交叉项的 XC 残差"""
    e_elec: float
    """完整电子能量"""
    e_nuc: float
    """核-核排斥能"""
    e_total: float
    """总能量 = e_elec + e_nuc"""

    @property
    def pi(self) -> float:
        """π 对电子能量的可加贡献 ``e_pi + e_cross``"""
        return self.e_pi + self.e_cross


def _to_pi(core) -> PiEnergy:  # type: ignore[no-untyped-def]
    return PiEnergy(
        e_pi=core.e_pi,
        e_sigma=core.e_sigma,
        e_cross=core.e_cross,
        e_cross_coul=core.e_cross_coul,
        e_cross_exch=core.e_cross_exch,
        e_cross_xc=core.e_cross_xc,
        e_elec=core.e_elec,
        e_nuc=core.e_nuc,
        e_total=core.e_total,
    )


class Calculator:
    """能量计算器

    基于分子轨道信息，对完整电子能量做 π/σ/πσ 分解。
    仅支持由高斯 log/out 读入、且为闭壳层、含轨道系数的分子。
    """

    def __init__(self, mole: Mole) -> None:
        self.mole = mole
        self.core = _core.props.moleprop.energy.Calculator(mole.core)  # type: ignore

    def pi_pocv_parts(self) -> PiEnergy:
        """POCV 方法的 π/σ/πσ 能量分解"""
        return _to_pi(self.core.pi_pocv_parts())

    def pi_mocv_parts(self) -> PiEnergy:
        """MOCV 方法的 π/σ/πσ 能量分解"""
        return _to_pi(self.core.pi_mocv_parts())

    def pi_pocv(self) -> float:
        """使用 POCV 方法计算的 π 对电子能量的贡献 ``E_π + E_πσ``"""
        return self.core.pi_pocv()

    def pi_mocv(self) -> float:
        """使用 MOCV 方法的 π 对电子能量的贡献 ``E_π + E_πσ``"""
        return self.core.pi_mocv()
