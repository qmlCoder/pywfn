"""
用来生成格点数据并保存cube文件
cube文件使用的是玻尔半径B
分子的坐标使用的是埃米A
A=B/1.889
B=A*1.889
计算的时候使用A，写文件的时候将坐标和格点转为B
"""

import numpy as np
from numpy.typing import NDArray

from pywfn import _core
from pywfn.base.mole import Mole


class CubWriter:
    """cube 格点文件写出器"""

    def __init__(self) -> None:
        """cube文件导出器"""
        self.core = _core.writer.CubWriter()  # type: ignore

    def set(
        self,
        title: str | None = None,  # 标题
        syms: list[str] | None = None,  # 元素符号
        xyzs: NDArray[np.float64] | None = None,  # 原子坐标
        obts: list[int] | None = None,  # 轨道序数
        pos0: list[float] | None = None,  # 起始坐标
        size: list[int] | None = None,  # 格点数量
        step: list[float] | None = None,  # 格点步长
        vals: NDArray[np.float64] | None = None,  # 格点数值[ngrid,nobt]
    ) -> None:
        """设置立方体文件的各个字段

        所有参数均为可选，只设置传入的字段；`obts` 与 `vals` 必须同时给出。

        Args:
            title: 标题字符串
            syms: 原子符号列表
            xyzs: 原子坐标列表
            obts: 原子核电荷数列表
            pos0: 立方体原点坐标 [x, y, z]
            size: 三个方向上的格点数 [nx, ny, nz]
            step: 三个方向上的步长 [dx, dy, dz]
            vals: 格点数据，形状为 (ngrid, nobt)
        """
        assert obts is not None and vals is not None, "obts和vals不能为None"
        assert len(obts) == vals.shape[0], "obts和vals的维度不匹配"
        self.core.set(title, syms, xyzs, obts, pos0, size, step, vals)

    @property
    def title(self) -> str:
        """标题字符串"""
        return self.core.title

    @title.setter
    def title(self, title: str) -> None:
        self.core.title = title

    @property
    def syms(self) -> list[str]:
        """原子符号列表"""
        return self.core.syms

    @syms.setter
    def syms(self, syms: list[str]) -> None:
        self.core.syms = syms

    @property
    def xyzs(self) -> list[list[float]]:
        """原子坐标列表"""
        return self.core.xyzs

    @xyzs.setter
    def xyzs(self, xyzs: list[list[float]]) -> None:
        self.core.xyzs = xyzs

    @property
    def obts(self) -> list[int]:
        """原子核电荷数列表"""
        return self.core.obts

    @obts.setter
    def obts(self, obts: list[int]) -> None:
        self.core.obts = obts

    @property
    def pos0(self) -> list[float]:
        """立方体原点坐标"""
        return self.core.pos0

    @pos0.setter
    def pos0(self, pos0: list[float]) -> None:
        self.core.pos0 = pos0

    @property
    def size(self) -> list[int]:
        """三个方向上的格点数"""
        return self.core.size

    @size.setter
    def size(self, size: list[int]) -> None:
        self.core.size = size

    @property
    def step(self) -> list[float]:
        """三个方向上的步长"""
        return self.core.step

    @step.setter
    def step(self, step: list[float]) -> None:
        self.core.step = step

    @property
    def vals(self) -> list[list[float]]:
        """三维格点数据"""
        return self.core.vals

    @vals.setter
    def vals(self, vals: NDArray[np.float64] | list[list[float]]) -> None:
        self.core.vals = vals

    def read_mole(self, mole: Mole) -> None:
        """从分子对象中读取数据"""
        self.core.read_mole(mole.core)

    def save(self, path: str) -> None:
        """保存cube文件"""
        self.core.save(path)
