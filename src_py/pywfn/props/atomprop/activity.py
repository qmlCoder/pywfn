"""
计算原子的反应活性，福井函数、parr函数等，一般需要多个分子才行
"""

import numpy as np

from pywfn import _core
from pywfn.base import Mole, Stm


class Calculator:
    def __init__(self, mole: Mole) -> None:
        self.mole = mole
        self.core = _core.props.atomprop.activity.Calculator(mole.core)  # type: ignore # 核心计算器

    def fukui(self, mole_n: Mole, mole_p: Mole, ctype: str = "mulliken") -> np.ndarray:
        """计算所有原子的福井函数

        Args:
            mole_n: 负电荷分子(N+1)
            mole_p: 正电荷分子(N-1)
            ctype: 电荷类型, 默认为 'mulliken'

        Returns:
            np.ndarray: 福井函数 [n,7] (N,N+1,N-1,f-,f+,f0,df)
        """
        vals = self.core.fukui(mole_n.core, mole_p.core, ctype)
        return np.array(vals)

    def fukui_pi_pocv(
        self, mole_n: Mole, mole_p: Mole, ctype: str = "mulliken"
    ) -> tuple[dict[int, list[float]], np.ndarray]:
        """基于pi电子数的福井函数"""
        dirs, vals = self.core.fukui_pi_pocv(mole_n.core, mole_p.core, ctype)
        return dirs, np.array(vals)

    def fukui_pi_mocv(
        self, mole_n: Mole, mole_p: Mole, ctype: str = "mulliken"
    ) -> tuple[dict[int, Stm], list[int], np.ndarray]:
        """基于pi电子数的福井函数"""
        stms, atos, vals = self.core.fukui_pi_mocv(mole_n.core, mole_p.core, ctype)
        stms = {idx: Stm(stm) for idx, stm in stms.items()}
        return stms, atos, np.array(vals)

    def fukui_dir(
        self,
        atm: int,
        dir: list[float],
        mole_n: Mole,
        mole_p: Mole,
        ctype: str = "mulliken",
    ) -> np.ndarray:
        """计算方向福井函数

        Args:
            atm: 要计算的原子索引
            dir: 原子的方向 [3]
            mole_n: 多一个电子的分子
            mole_p: 少一个电子的分子
            ctype: 电荷类型, 默认为 mulliken

        Returns:
            np.ndarray: 指定方向的福井函数 [n,6] (N,N+1,N-1,f-,f+,f0,df)
        """
        vals = self.core.fukui_dir(atm, dir, mole_n.core, mole_p.core, ctype)
        return np.array(vals)

    def freev(self, otype: str = "mayer") -> np.ndarray:
        """计算所有原子的自由价

        Args:
            otype: 键级类型, 默认为 "mayer"

        Returns:
            np.ndarray: 形状为 [n, 3]，每行为 (标准价, 总键级, 自由价)
        """
        return np.array(self.core.freev(otype))

    def freev_pi(
        self, otype: str = "mayer"
    ) -> tuple[dict[int, list[float]], np.ndarray]:
        """使用pi键级计算的自由价

        Args:
            otype: 键级类型, 默认为 "mayer"

        Returns:
            tuple: (方向字典, 自由价值)
        """
        dirs, vals = self.core.freev_pi(otype)
        return dirs, np.array(vals)

    def freev_dir(
        self,
        atm: int,
        dir: list[float],
        nebs: list[int],
        otype: str = "mayer",
        mtype: str = "pocv",
    ) -> tuple[float, np.ndarray, float]:
        """计算指定方向上的自由价

        Args:
            atm: 原子索引
            dir: 方向向量 [3]
            nebs: 相邻原子索引列表
            otype: 键级类型, 默认为 "mayer"
            mtype: 分解方法, 可选 "pocv" 或 "mocv", 默认为 "pocv"

        Returns:
            tuple: (标准价, 各相邻原子在此方向的键级, 方向自由价)
        """
        stand, orders, value = self.core.freev_dir(atm, dir, nebs, otype, mtype)
        return stand, np.array(orders), value
