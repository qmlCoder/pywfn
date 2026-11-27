from dataclasses import dataclass

from pywfn import _core
from pywfn.base.mole import Mole


@dataclass
class BondEnergy:
    """键能量的可加分解结果（闭壳层，杂化泛函）

    满足 ``e_bond + e_rest + e_cross == e_elec``。
    """

    e_bond: float
    """键自能量"""
    e_rest: float
    """其余部分能量"""
    e_cross: float
    """键-其余交叉项"""
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
    def bond(self) -> float:
        """键对电子能量的可加贡献 ``e_bond + e_cross``"""
        return self.e_bond + self.e_cross


def _to_bond(core) -> BondEnergy:  # type: ignore[no-untyped-def]
    return BondEnergy(
        e_bond=core.e_bond,
        e_rest=core.e_rest,
        e_cross=core.e_cross,
        e_cross_coul=core.e_cross_coul,
        e_cross_exch=core.e_cross_exch,
        e_cross_xc=core.e_cross_xc,
        e_elec=core.e_elec,
        e_nuc=core.e_nuc,
        e_total=core.e_total,
    )


class Calculator:
    """键能量计算器

    把完整电子能量按"键 / 其余 / 交叉"做可加分解。
    仅支持由高斯 log/out 读入、且为闭壳层、含轨道系数的分子。
    """

    def __init__(self, mole: Mole) -> None:
        self.mole = mole
        self.caler = _core.props.bondprop.energy.Calculator(mole.core)  # type: ignore

    def bond_energy(self, atm1: int, atm2: int) -> float:
        """键（两原子）对电子能量的可加贡献 ``E_bond + E_cross``"""
        return self.caler.bond_energy(atm1, atm2)

    def bond_parts(self, atm1: int, atm2: int) -> BondEnergy:
        """对指定原子对的完整电子能量做 bond / rest / cross 分解"""
        return _to_bond(self.caler.bond_parts(atm1, atm2))
