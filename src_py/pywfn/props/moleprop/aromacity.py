"""
计算分子的芳香性
使用pi键级的标准差表示
"""

import numpy as np

from pywfn import _core
from pywfn.base.mole import Mole


class Calculator:
    def __init__(self, mole: Mole):
        self.mole = mole
        self.ratio = 0.5
        self.caler = _core.props.moleprop.aromacity.Calculator(mole.core)  # type: ignore

    def PISD(self, rings: list[list[int]] | None = None) -> np.ndarray:
        """使用pi键级标准差计算芳香性

        Args:
            rings: 环的原子索引列表, 如 [[1,2,3,4,5,6], ...]
                    如果为None则自动检测所有环

        Returns:
            np.ndarray: 每个环的芳香性指数
        """
        return self.caler.PISD(rings)
