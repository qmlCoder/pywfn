"""
定义一个平面的方式有很多种，但是本质上对应的格点都是一样的，因此可以用统一的类来实现
"""

import numpy as np

from pywfn import _core


class LineGrid:
    def __init__(self, core: _core.props.gridprop.LineGrid):  # type: ignore
        """线网格"""
        self.core = core

    def set(self, p0: np.ndarray, p1: np.ndarray, step: float) -> None:
        """设置线网格参数

        Args:
            p0: 起点坐标
            p1: 终点坐标
            step: 步长
        """
        self.core.set(p0, p1, step)

    def get(self) -> np.ndarray:
        """获取格点坐标"""
        return self.core.get()


class RectGrid:
    def __init__(self, core: _core.props.gridprop.RectGrid):  # type: ignore
        """矩形网格"""
        self.core = core

    def set_v1(
        self,
        center: np.ndarray,
        normal: np.ndarray,
        vector: np.ndarray,
        size: float,
        step: float,
    ) -> None:
        """通过中心、法向量和方向向量设置矩形网格

        Args:
            center: 中心点坐标
            normal: 法向量
            vector: 方向向量
            size: 网格大小
            step: 步长
        """
        self.core.set_v1(center, normal, vector, size, step)

    def set_v2(
        self, p0: np.ndarray, p1: np.ndarray, p2: np.ndarray, step: float
    ) -> None:
        """根据三个原子坐标生成矩形网格"""
        self.core.set_v2(p0, p1, p2, step)

    def get(self) -> tuple[list[int], np.ndarray]:
        """获取格点坐标"""
        return self.core.get()


class CubeGrid:
    def __init__(self):
        """立方体网格"""
        self.core = _core.props.gridprop.CubeGrid()  # type: ignore

    def set_v1(
        self, p0: list[float], p1: list[float], step: float, bord: float = 0.0
    ) -> None:
        """设置立方体网格参数

        Args:
            p0: 起始点坐标
            p1: 终止点坐标
            step: 步长
            bord: 边界扩展
        """
        self.core.set_v1(p0, p1, step, bord)

    def get(self) -> tuple[list[int], np.ndarray]:
        """获取格点坐标和数量"""
        return self.core.get()
