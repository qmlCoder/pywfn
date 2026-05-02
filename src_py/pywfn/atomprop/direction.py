"""
计算原子的各种方向
- 原子法向量
- p轨道方向
- 原子周围一圈方向
"""
from pywfn.base.mole import Mole
from pywfn import _core
from pywfn.base import Stm
import numpy as np

maxWeaves = {}  # 记录已经计算过的原子最大值方向


class Calculator:  # type: ignore
    def __init__(self, mole: Mole) -> None:
        self.core = _core.atomprop.direction.Calculator(mole.core) # type: ignore

    def normal_vector(self, atm: int):
        """计算原子的法向量"""
        return self.core.normal_vector(atm)

    def get_normal(self, atm: int):
        """计算原子的法向量"""
        return self.core.get_normal(atm)

    def LCS(self, atm: int, neb: int | None = None):
        """计算原子的局部坐标系"""
        return Stm(self.core.LCS(atm, neb))

    def reactions(self, atm: int):
        """计算原子可能的反应方向"""
        res=self.core.reactions(atm)
        return np.array(res)
