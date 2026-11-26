"""
三维向量对象
封装 rswfn 的 Vec3，提供常用的向量运算

- 创建：new / from_arr / from_vec / zero / rand
- 分量：x / y / z
- 运算：加减、数乘、数除、取负
- 几何：长度、归一化、叉积、点积、夹角、投影、旋转、点到平面的距离
"""

import numpy as np

from pywfn import _core


class Vec3:
    """三维向量"""

    def __init__(self, _core: _core.base.Vec3) -> None:  # type: ignore
        self.core = _core

    @staticmethod
    def new(x: float, y: float, z: float) -> "Vec3":
        """创建一个新的三维向量

        Args:
            x: x 分量
            y: y 分量
            z: z 分量
        """
        return Vec3(_core.base.Vec3(x, y, z))  # type: ignore

    @staticmethod
    def from_arr(vec: list[float]) -> "Vec3":
        """从长度为3的数组创建三维向量

        Args:
            vec: [x, y, z] 数组
        """
        return Vec3(_core.base.Vec3.from_arr(vec))  # type: ignore

    @staticmethod
    def from_vec(vec: list[float]) -> "Vec3":
        """从列表中创建三维向量

        Args:
            vec: 至少包含3个元素的列表
        """
        return Vec3(_core.base.Vec3.from_vec(vec))  # type: ignore

    @staticmethod
    def zero() -> "Vec3":
        """返回零向量 (0, 0, 0)"""
        return Vec3(_core.base.Vec3.zero())  # type: ignore

    @staticmethod
    def rand() -> "Vec3":
        """生成一个随机单位向量"""
        return Vec3(_core.base.Vec3.rand())  # type: ignore

    @property
    def x(self) -> float:
        """x 分量"""
        return self.core.x

    @property
    def y(self) -> float:
        """y 分量"""
        return self.core.y

    @property
    def z(self) -> float:
        """z 分量"""
        return self.core.z

    def into_arr(self) -> list[float]:
        """将向量转换为长度为3的数组"""
        return self.core.into_arr()

    def into_vec(self) -> list[float]:
        """将向量转换为列表"""
        return self.core.into_vec()

    def normalize(self) -> "Vec3":
        """返回归一化后的单位向量"""
        return Vec3(self.core.normalize())

    def length(self) -> float:
        """返回向量的长度（模）"""
        return self.core.length()

    def cross(self, vec3: "Vec3") -> "Vec3":
        """计算与另一个向量的叉积

        Args:
            vec3: 另一个三维向量
        """
        return Vec3(self.core.cross(vec3.core))

    def dot(self, vec3: "Vec3") -> float:
        """计算与另一个向量的点积

        Args:
            vec3: 另一个三维向量
        """
        return self.core.dot(vec3.core)

    def angle(self, vec3: "Vec3") -> float:
        """计算与另一个向量之间的夹角，范围 0-PI

        Args:
            vec3: 另一个三维向量
        """
        return self.core.angle(vec3.core)

    def project_onto(self, vec3: "Vec3") -> "Vec3":
        """计算向量在另一个向量方向上的投影

        Args:
            vec3: 目标投影方向向量
        """
        return Vec3(self.core.project_onto(vec3.core))

    def rotate_arrow(self, axis: "Vec3", angle: float) -> "Vec3":
        """绕指定轴旋转，self 代表箭头的朝向

        Args:
            axis: 旋转轴
            angle: 旋转角度
        """
        return Vec3(self.core.rotate_arrow(axis.core, angle))

    def rotate_point(self, cent: "Vec3", axis: "Vec3", angle: float) -> "Vec3":
        """绕指定中心点和轴旋转，self 代表点

        Args:
            cent: 旋转中心
            axis: 旋转轴
            angle: 旋转角度
        """
        return Vec3(self.core.rotate_point(cent.core, axis.core, angle))

    def dist_to_plan(self, plan: list[float]) -> float:
        """计算点到平面的距离，平面由 ax + by + cz + d = 0 定义

        Args:
            plan: 平面参数 [a, b, c, d]
        """
        return self.core.dist_to_plan(plan)

    def __add__(self, vec3: "Vec3") -> "Vec3":
        return Vec3(self.core + vec3.core)

    def __radd__(self, vec3: "Vec3") -> "Vec3":
        return Vec3(vec3.core + self.core)

    def __sub__(self, vec3: "Vec3") -> "Vec3":
        return Vec3(self.core - vec3.core)

    def __rsub__(self, vec3: "Vec3") -> "Vec3":
        return Vec3(vec3.core - self.core)

    def __mul__(self, scalar: float) -> "Vec3":
        return Vec3(self.core * scalar)

    def __rmul__(self, scalar: float) -> "Vec3":
        return Vec3(scalar * self.core)

    def __truediv__(self, scalar: float) -> "Vec3":
        return Vec3(self.core / scalar)

    def __neg__(self) -> "Vec3":
        return Vec3(-self.core)

    def __array__(self, dtype=None, copy=None) -> np.ndarray:
        return np.array(self.core.into_arr(), dtype=dtype)

    def __repr__(self) -> str:
        return self.core.__repr__()
