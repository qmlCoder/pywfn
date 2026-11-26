"""
计算DFT格点

镜像自 Rust 模块 `src_rs/gridprop/dftgrid.rs` 的结构，不包含实际代码逻辑。
"""

from typing import Tuple

from numpy import ndarray

from pywfn import _core
from pywfn.base.mole import Mole


class Calculator:
    def __init__(self, mole: Mole) -> None:
        self.mole = mole
        self.core = _core.props.gridprop.dftgrid.Calculator(mole.core)  # type: ignore

    def atm_grids(
        self, iatm: int, nrad: int, fsph: int
    ) -> tuple[list[list[float]], list[float]]:
        """计算指定原子的 DFT 格点坐标及权重。

        Args:
            iatm: 原子索引（从 0 开始）
            nrad: 径向格点数量
            fsph: 球面格点函数编号，控制球面积分精度

        Returns:
            (格点坐标, 格点权重) 二元组
        """
        return self.core.atm_grids(iatm, nrad, fsph)

    def a2m_weits(
        self,
        iatm: int,
        atm_grids: list[list[float]],
        atm_weits: list[float],
    ) -> list[float]:
        """将原子格点权重转换为分子格点权重（Becke 权重）。

        Args:
            iatm: 原子索引（从 0 开始）
            atm_grids: 该原子的格点坐标列表
            atm_weits: 该原子的初始格点权重

        Returns:
            Becke 权重列表，用于分子数值积分
        """
        return self.core.a2m_weits(iatm, atm_grids, atm_weits)

    def mol_grids(self, nrad: int, nsph: int) -> Tuple[ndarray, ndarray]:
        """计算整个分子的 DFT 格点及权重。

        遍历所有原子生成径向-球面格点，计算 Becke 权重后合并为分子格点。

        Args:
            nrad: 径向格点数量
            nsph: 球面格点数量 可选值为: 6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434, 590

        Returns:
            (分子格点坐标, 分子格点权重) 二元组
            坐标数组形状为 (N, 3)，权重数组形状为 (N,)
        """
        return self.core.mol_grids(nrad, nsph)

    def frag_grids(self, frag: list[int]) -> tuple[list[list[float]], list[float]]:
        """计算分子片段（指定原子集合）的 DFT 格点及权重。

        Args:
            frag: 片段包含的原子索引列表

        Returns:
            (片段格点坐标, 片段格点权重) 二元组
        """
        return self.core.frag_grids(frag)


def rad_grids(R: float, nrad: int) -> tuple[list[float], list[float]]:
    """计算 DFT 原子径向格点及其权重。

    Args:
        R: 原子半径
        nrad: 径向格点数量

    Returns:
        (径向格点, 径向权重) 二元组
    """
    return _core.props.gridprop.dftgrid.rad_grids(R, nrad)  # type: ignore


def sph_grids(
    nsph: int,
) -> tuple[ndarray, ndarray]:
    """计算 DFT 球面格点及其权重。

    Args:
        nsph: 球面格点数量 可选值为 6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434, 590

    Returns:
        (球面格点, 球面权重) 二元组
    """
    return _core.props.gridprop.dftgrid.sph_grids(nsph)  # type: ignore
