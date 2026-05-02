"""
定义一个平面的方式有很多种，但是本质上对应的格点都是一样的，因此可以用统一的类来实现
"""

import numpy as np
from numpy.typing import NDArray
from pywfn.base.mole import Mole
from pywfn import _core


class LineGrid:
    def __init__(self, core: _core.gridprop.LineGrid):  # type: ignore
        self.core = core

    def set(self, p0: np.ndarray, p1: np.ndarray, step: float):
        self.core.set(p0, p1, step)

    def get(self):
        return self.core.get()


class RectGrid:
    def __init__(self, core: _core.gridprop.RectGrid):  # type: ignore
        self.core = core

    def set_v1(self, center: np.ndarray, normal: np.ndarray, vector: np.ndarray, size: float, step: float):
        self.core.set_v1(center, normal, vector, size, step)

    def set_v2(self, p0: np.ndarray, p1: np.ndarray, p2: np.ndarray, step: float):  # 根据三个原子的坐标生成网格
        self.core.set_v2(p0, p1, p2, step)

    def get(self):
        return self.core.get()


class CubeGrid:
    def __init__(self):  # type: ignore
        self.core = _core.gridprop.CubeGrid() # type: ignore

    def set_v1(self, p0: list[float], p1: list[float], step: float, bord: float = 0.0):
        self.core.set_v1(p0, p1, step, bord)

    def get(self):
        return self.core.get()
