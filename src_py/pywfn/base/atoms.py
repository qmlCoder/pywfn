import numpy as np

import pywfn
from pywfn import _core
from pywfn.base.vect import Vec3


class Stm:
    """原子局部坐标系"""

    def __init__(self, _core: _core.base.Stm) -> None:  # type: ignore
        self.core = _core

    @staticmethod
    def new(ex: list[float], ey: list[float], ez: list[float]):
        core = _core.base.Stm(ex, ey, ez)  # type: ignore
        return Stm(core)

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
        return f"{self.core.__repr__()}"


class Atom:
    """
    原子对象
    """

    def __init__(self, _core: _core.base.Atom) -> None:  # type: ignore
        self.core = _core

    @property
    def rad(self) -> float:
        """原子半径"""
        return self.core.rad

    @property
    def sym(self) -> str:
        """元素符号"""
        return self.core.sym


class Atoms:  # type: ignore
    def __init__(self, core: _core.base.Atoms, mole: "pywfn.base.Mole") -> None:  # type: ignore
        self.core = core
        self.mole = mole

    def add(self, atm: int, xyz: list[float], crg: float) -> None:
        """添加一个原子"""
        self.core.add(atm, xyz, crg)

    def get(self, idx: int) -> Atom:
        """获取指定索引的原子"""
        return Atom(self.core.get(idx))

    def len(self) -> int:
        """原子数量"""
        return self.core.len()

    def vec(self, ai: int, aj: int) -> Vec3:
        """计算 ai 指向 aj 的向量"""
        return Vec3(self.core.vec(ai, aj))

    @property
    def neles(self):
        return self.core.neles()

    def dist(self, i: int, j: int) -> float:
        """计算两个原子之间的距离"""
        return self.core.dist(i, j)

    @property
    def syms(self) -> list[str]:
        """所有原子的元素符号"""
        return self.core.syms()

    @property
    def xyzs(self) -> np.ndarray:
        """所有原子的坐标"""
        return np.array(self.core.xyzs())

    def nebs(self, atm: int) -> list[int]:
        """获取指定原子的相邻原子索引"""
        return self.core.nebs(atm)

    @property
    def nucs(self) -> list[int]:
        """所有原子的核电荷数"""
        return self.core.nucs()

    @property
    def crgs(self) -> list[float]:
        """所有原子的电荷"""
        return self.core.crgs()

    def is_line(self, atm: int) -> bool:
        """判断原子是否在一条直线上"""
        return self.core.is_line(atm)

    def is_bond(self, i: int, j: int) -> bool:
        """判断两个原子之间是否成键"""
        return self.core.is_bond(i, j)

    def __repr__(self) -> str:
        return self.core.__repr__()

    def bond_len(self,ai:int,aj:int)->float:
        """计算键长"""
        return self.core.bond_len(ai,aj)

    def bond_ang(self,ai:int,aj:int,ak:int)->float:
        """计算键角"""
        return self.core.bond_ang(ai,aj,ak)

    def face_ang(self,ai:int,aj:int,ak:int,al:int)->float:
        """计算二面角"""
        return self.core.face_ang(ai,aj,ak,al)
