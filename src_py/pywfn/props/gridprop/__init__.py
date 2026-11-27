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

    def set_v3(
        self, p00: np.ndarray, p01: np.ndarray, p11: np.ndarray, size: list[int]
    ) -> None:
        """由三个顶点和格点数量确定矩形网格

        Args:
            p00: 第一个顶点
            p01: 第二个顶点
            p11: 第三个顶点
            size: 两个方向的格点数量
        """
        self.core.set_v3(p00, p01, p11, size)

    def paras(
        self, p00: np.ndarray, p01: np.ndarray, p11: np.ndarray
    ) -> list[float]:
        """计算矩形网格的平面参数

        Args:
            p00: 第一个顶点
            p01: 第二个顶点
            p11: 第三个顶点

        Returns:
            list[float]: [a, b, c, d] 四个参数
        """
        return self.core.paras(p00, p01, p11)

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

class GridVals:
    """格点上函数值及其各级导数"""

    def __init__(
        self,
        val0: float = 0.0,
        val1: list[float] | None = None,
        val2: list[list[float]] | None = None,
    ) -> None:
        self.core = _core.props.gridprop.GridVals(val0, val1, val2)  # type: ignore

    @staticmethod
    def non() -> "GridVals":
        """返回一个全零的格点函数值对象"""
        gv = GridVals.__new__(GridVals)
        gv.core = _core.props.gridprop.GridVals.non()  # type: ignore
        return gv

    def set(self, val0: float, val1: list[float], val2: list[list[float]]) -> None:
        """设置函数值及其各级导数

        Args:
            val0: 函数值
            val1: 梯度 [dx, dy, dz]
            val2: Hessian 矩阵 3x3
        """
        self.core.set(val0, val1, val2)

    @property
    def val0(self) -> float:
        """函数值"""
        return self.core.val0

    @property
    def val1(self) -> list[float]:
        """梯度"""
        return self.core.val1

    @property
    def val2(self) -> list[list[float]]:
        """Hessian 矩阵"""
        return self.core.val2

    def grad_norm(self) -> float:
        """梯度的模"""
        return self.core.grad_norm()

    def grad_norm2(self) -> float:
        """梯度模的平方"""
        return self.core.grad_norm2()

    def val2s(self) -> list[float]:
        """Hessian 矩阵展平为长度 9 的列表"""
        return self.core.val2s()