import numpy as np

import pywfn
from pywfn import _core


class Stm:
    """原子局部坐标系"""

    def __init__(self, _core: _core.base.Stm) -> None:  # type: ignore
        self.core = _core

    def set(self, ex: list[float], ey: list[float], ez: list[float]) -> None:
        self.core.set(ex, ey, ez)

    @property
    def ex(self) -> list[float]:
        return self.core.ex

    @property
    def ey(self) -> list[float]:
        return self.core.ey

    @property
    def ez(self) -> list[float]:
        return self.core.ez

    @property
    def mat(self) -> np.ndarray:
        return self.core.mat

    def __repr__(self) -> str:
        return f"pywfn: {self.core.__repr__()}"


class Atom:
    """
    原子对象
    """

    def __init__(self, _core: _core.base.Atom) -> None:  # type: ignore
        self.core = _core

    @property
    def rad(self):
        return self.core.rad

    @property
    def sym(self):
        return self.core.sym


class Atoms:  # type: ignore
    def __init__(self, core: _core.base.Atoms, mole: "pywfn.base.Mole") -> None:  # type: ignore
        self.core = core
        self.mole = mole

    def len(self):
        return self.core.len()

    def dist(self, i: int, j: int):
        return self.core.dist(i, j)

    @property
    def syms(self):
        return self.core.syms()

    @property
    def xyzs(self):
        return self.core.xyzs()

    def nebs(self, atm: int) -> list[int]:
        return self.core.nebs(atm)

    @property
    def nucs(self):
        return self.core.nucs()

    @property
    def crgs(self):
        return self.core.crgs()

    def is_line(self, atm: int) -> bool:
        return self.core.is_line(atm)

    def is_bond(self, i: int, j: int) -> bool:
        return self.core.is_bond(i, j)

    def __repr__(self) -> str:
        return self.core.__repr__()
