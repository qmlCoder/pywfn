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

    def check(self) -> bool:
        """判断变换矩阵是否有效"""
        return self.core.check()

    @staticmethod
    def blank() -> "Stm":
        """创建一个空的变换矩阵"""
        return Stm(_core.base.Stm.blank())  # type: ignore

    def vals(self) -> list[float]:
        """把矩阵展平为长度 9 的列表"""
        return self.core.vals()

    @staticmethod
    def from_mat(mat: np.ndarray) -> "Stm":
        """从 3x3 矩阵创建变换矩阵

        Args:
            mat: 3x3 矩阵
        """
        return Stm(_core.base.Stm.from_mat(mat))  # type: ignore

    def relative(self, stmb: "Stm") -> "Stm":
        """计算相对另一个变换矩阵的旋转

        Args:
            stmb: 另一个变换矩阵
        """
        return Stm(self.core.relative(stmb.core))

    def dirs(self) -> list[list[float]]:
        """获取变换矩阵的三个行向量"""
        return self.core.dirs()

    def copy(self, stm: "Stm") -> None:
        """从另一个变换矩阵复制内容

        Args:
            stm: 源变换矩阵
        """
        self.core.copy(stm.core)

    @staticmethod
    def rand_ex(ex: list[float]) -> "Stm":
        """由给定 x 轴构造随机变换矩阵

        Args:
            ex: x 轴方向
        """
        return Stm(_core.base.Stm.rand_ex(ex))  # type: ignore

    @staticmethod
    def rand_ez(ez: list[float]) -> "Stm":
        """由给定 z 轴构造随机变换矩阵

        Args:
            ez: z 轴方向
        """
        return Stm(_core.base.Stm.rand_ez(ez))  # type: ignore

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

    def vec(self) -> Vec3:
        """由坐标原点指向该原子的向量"""
        return Vec3(self.core.vec())

    def irow(self) -> int:
        """原子在周期表中的行号"""
        return self.core.irow()


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

    def dmat(self) -> np.ndarray:
        """所有原子间的距离矩阵"""
        return self.core.dmat()

    def nuc_eng(self) -> float:
        """核间排斥能"""
        return self.core.nuc_eng()

    def xyzs_ang(self) -> np.ndarray:
        """所有原子的坐标（单位为埃）"""
        return np.array(self.core.xyzs_ang())

    def rads(self) -> list[float]:
        """所有原子的半径列表"""
        return self.core.rads()

    def stms(self) -> list["Stm | None"]:
        """所有原子的局域坐标系"""
        return [Stm(stm) if stm is not None else None for stm in self.core.stms()]

    def idxs(self) -> list[int]:
        """所有原子的索引列表"""
        return self.core.idxs()

    def ang(self, ai: int, aj: int, ak: int) -> float:
        """计算三个原子形成的夹角

        Args:
            ai: 第一个原子索引
            aj: 顶点原子索引
            ak: 第三个原子索引
        """
        return self.core.ang(ai, aj, ak)

    def bonds(self) -> list[list[int]]:
        """分子中的成键原子对"""
        return self.core.bonds()

    def set_bonds(self, bonds: list[list[int]]) -> None:
        """设置分子中的成键原子对

        Args:
            bonds: 成键原子对列表
        """
        self.core.set_bonds(bonds)

    def cent(self, atms: list[int]) -> Vec3:
        """计算指定原子集合的几何中心

        Args:
            atms: 原子索引列表
        """
        return Vec3(self.core.cent(atms))

    def equal(self, other: "Atoms") -> bool:
        """判断两个原子集合是否完全相同

        Args:
            other: 另一个原子集合
        """
        return self.core.equal(other.core)
