"""
定义分子结构编辑器
可以用来生成刚性扫描的结构
"""

import numpy as np
from numpy.typing import NDArray

from pywfn import _core
from pywfn.base import Mole


def rotate_points(
    points: NDArray[np.float64],
    cent: NDArray[np.float64],
    axis: NDArray[np.float64],
    angle: float,
) -> NDArray[np.float64]:
    """绕过指定中心点和轴的旋转所有点

    Args:
        points: 待旋转的点坐标，形状为 (n, 3)
        cent: 旋转中心
        axis: 旋转轴方向
        angle: 旋转角度（弧度）

    Returns:
        NDArray[np.float64]: 旋转后的点坐标，形状为 (n, 3)
    """
    return np.array(
        _core.editor.rotate_points(points, cent, axis, angle)  # type: ignore
    )


class Editor:
    """分子结构编辑器"""

    def __init__(self, mole: Mole) -> None:
        self.mole = mole
        self.core = _core.editor.Editor(mole.core)  # type: ignore

    def get_group(self, atm: int, excludes: list[int] | None = None) -> list[int]:
        """获取与指定原子直接相连的原子分组

        Args:
            atm: 目标原子索引
            excludes: 需要排除的原子索引

        Returns:
            list[int]: 相连原子的索引列表
        """
        return self.core.get_group(atm, excludes or [])

    def set_bl_by(self, a1: int, a2: int, value: float, stype: int = 0) -> Mole:
        """按键长设置分子结构

        Args:
            a1: 第一个原子索引
            a2: 第二个原子索引
            value: 目标键长
            stype: 结构类型

        Returns:
            Mole: 修改后的新分子对象
        """
        return Mole(self.core.set_bl_by(a1, a2, value, stype))

    def set_ba_by(
        self, a1: int, a2: int, a3: int, value: float, btype: int = 0
    ) -> Mole:
        """按键角设置分子结构

        Args:
            a1: 第一个原子索引
            a2: 第二个原子索引
            a3: 第三个原子索引
            value: 目标键角（弧度）
            btype: 结构类型

        Returns:
            Mole: 修改后的新分子对象
        """
        return Mole(self.core.set_ba_by(a1, a2, a3, value, btype))

    def set_da_by(self, a1: int, a2: int, value: float, rtype: int = 0) -> Mole:
        """按二面角设置分子结构

        Args:
            a1: 第一个原子索引
            a2: 第二个原子索引
            value: 目标二面角（弧度）
            rtype: 结构类型

        Returns:
            Mole: 修改后的新分子对象
        """
        return Mole(self.core.set_da_by(a1, a2, value, rtype))
