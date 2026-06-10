from typing import Literal

from pywfn import _core
from pywfn.base import Mole, Stm

Chrgs = Literal["mulliken", "lowdin", "space", "hirshfeld"]


class Calculator:
    def __init__(self, mole: "Mole", atms: list[int] | None = None):
        self.mole = mole
        self.caler = _core.atomprop.charge.Calculator(mole.core, atms)  # type: ignore # 核心计算器

    def spin(self, ctype: str = "mulliken"):
        """计算所有原子的自旋"""
        return self.caler.spin(ctype)

    def mulliken(self):
        """
        计算mulliken电荷
        """
        return self.caler.mulliken()

    def lowdin(self):
        """
        计算每个原子的lowdin电荷
        """
        return self.caler.lowdin()

    def hirshfeld(self):
        """计算原子的Hirshfeld电荷"""
        return self.caler.hirshfeld()

    def pocv(
        self,
        dirs: dict[int, list[float]],
        keep_other_atm: bool,
        keep_other_sym: bool,
        ctype: str,
    ):
        """计算投影分子轨道的电子数"""
        return self.caler.pocv(dirs, keep_other_atm, keep_other_sym, ctype)  # type: ignore

    def mocv(self, stms: dict[int, Stm], atos: list[int], ctype: str = "mulliken"):
        """
        计算mocv电荷
        """
        stms = {atm: stm.core for atm, stm in stms.items()}
        return self.caler.mocv(stms, atos, ctype)

    def pi_pocv(self, ctype: str = "mulliken"):
        """
        计算π电子数
        每个原子的方向为`法向量`的`方向电子数`即为`π电子数`
        """
        return self.caler.pi_pocv(ctype)

    def pi_mocv(
        self, ctype="mulliken"
    ):  # 使用轨道分解方法计算pi电子分布，可以包含D轨道
        return self.caler.pi_mocv(ctype)


def mulliken(dmat, smat, uplo):
    return _core.atomprop.charge.mulliken(dmat, smat, uplo)  # type: ignore


def lowdin(dmat, smat, uplo):
    return _core.atomprop.charge.lowdin(dmat, smat, uplo)  # type: ignore
