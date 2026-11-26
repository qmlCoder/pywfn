"""
计算与键相关的方向
"""

import numpy as np

from pywfn import _core
from pywfn.base.mole import Mole


class Calculator:
    def __init__(self, mole: Mole) -> None:
        self.core = _core.props.bondprop.direction.Calculator(mole.core)  # type: ignore

    def verts(self, ai: int, aj: int) -> np.ndarray:
        """计算垂直于某根键的方向

        Args:
            ai: 第一个原子索引
            aj: 第二个原子索引

        Returns:
            np.ndarray: 垂直于键的方向向量
        """
        dirs = self.core.verts(ai, aj)
        return np.array(dirs)
