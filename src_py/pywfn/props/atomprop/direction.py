"""
计算原子的各种方向
- 原子法向量
- p轨道方向
- 原子周围一圈方向
"""

import numpy as np

from pywfn import _core
from pywfn.base import Stm
from pywfn.base.mole import Mole

maxWeaves = {}  # 记录已经计算过的原子最大值方向


class Calculator:  # type: ignore
    def __init__(self, mole: Mole) -> None:
        self.core = _core.props.atomprop.direction.Calculator(mole.core)  # type: ignore

    def normal_vector(self, atm: int) -> np.ndarray:
        """计算原子的法向量"""
        return self.core.normal_vector(atm)

    def get_normal(self, atm: int) -> np.ndarray:
        """计算原子的法向量"""
        return self.core.get_normal(atm)

    def LCS(self, atm: int, neb: int | None = None) -> Stm | None:
        """计算原子的局部坐标系"""
        stm = self.core.LCS(atm, neb)
        return Stm(stm) if stm is not None else None

    def reactions(self, atm: int) -> np.ndarray:
        """计算原子可能的反应方向"""
        res: list[list[float]] = self.core.reactions(atm)
        return np.array(res)

    def dmat_fix(self) -> np.ndarray:
        """返回修正后的密度矩阵"""
        return self.core.dmat_fix()

    def RLCS(self, atm: int) -> Stm:
        """使用随机局域坐标系（RLCS）

        Args:
            atm: 原子索引
        """
        return Stm(self.core.RLCS(atm))

    def MLCS(self, atm: int) -> list[Stm]:
        """使用多个局域坐标系（MLCS）

        Args:
            atm: 原子索引
        """
        return [Stm(stm) for stm in self.core.MLCS(atm)]

    def gt_bond_angle(self, atm: int) -> np.ndarray:
        """获取指定原子的键角方向

        Args:
            atm: 原子索引
        """
        return np.array(self.core.gt_bond_angle(atm))

    def linear_circle(self, atm: int) -> np.ndarray:
        """获取线性分子的圆周方向

        Args:
            atm: 原子索引
        """
        return np.array(self.core.linear_circle(atm))

    def min_exclusion(self, atm: int) -> np.ndarray:
        """获取排除其他原子后的最小方向

        Args:
            atm: 原子索引
        """
        return self.core.min_exclusion(atm)
