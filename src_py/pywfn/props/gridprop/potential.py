import numpy as np

from pywfn import _core
from pywfn.base.mole import Mole


class Calculator:
    def __init__(self, mole: Mole) -> None:
        self.mole = mole
        self.core = _core.props.gridprop.potential.Calculator(mole.core)  # type: ignore

    def nuc_potential(self, grids: np.ndarray) -> list[float]:
        """计算原子核势能（核-电子吸引势）在给定格点上的值

        Args:
            grids: 三维坐标点列表 [[x, y, z], ...]

        Returns:
            每个格点上的核势能值
        """
        return self.core.nuc_potential(grids)

    def ele_potential(self, grids: np.ndarray) -> list[float]:
        """计算电子势能（电子-电子排斥势）在给定格点上的值

        Args:
            grids: 三维坐标点列表 [[x, y, z], ...]

        Returns:
            每个格点上的电子势能值
        """
        return self.core.ele_potential(grids)
