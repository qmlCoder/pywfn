"""
计算分子空间电子密度
整个分子只有网格点，但是没有权重
对于每个原子，根据自身的坐标及权重插值到整个空间的坐标及权重
？ 波函数值是否也要插值呢？

空间中一个格点的权重是否应该对所有分子都一致？
"""

import numpy as np

from pywfn import _core
from pywfn.base.mole import Mole


class Calculator:
    def __init__(self, mole: Mole) -> None:
        self.caler = _core.props.gridprop.density.Calculator(mole.core)  # type: ignore # 核心计算器

    def mol_rho_cm(
        self, grids: np.ndarray, level: int
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """以分子轨道的方式计算电子密度，可计算每个分子轨道的电子密度

        Args:
            grids: 三维坐标点列表 [[x, y, z], ...]
            level: 计算级别（0=仅密度，1=密度+梯度，2=密度+梯度+Hessian）

        Returns:
            (密度值, 梯度, Hessian) 三元组
        """
        return self.caler.mol_rho_cm(grids, level)

    def mol_rho_dm(
        self, grids: np.ndarray, level: int
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """以密度矩阵的方式计算电子密度（适用于开壳层）

        Args:
            grids: 三维坐标点列表 [[x, y, z], ...]
            level: 计算级别（0=仅密度，1=密度+梯度，2=密度+梯度+Hessian）

        Returns:
            (密度值, 梯度, Hessian) 三元组
        """
        return self.caler.mol_rho_dm(grids, level)

    def pi_rho_pocv(
        self, grids: np.ndarray, level: int
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """计算 pi 电子密度（POCV）

        Args:
            grids: 三维坐标点列表 [[x, y, z], ...]
            level: 计算级别（0=仅密度，1=密度+梯度，2=密度+梯度+Hessian）

        Returns:
            (密度值, 梯度, Hessian) 三元组
        """
        return self.caler.pi_rho_pocv(grids, level)

    def pi_rho_mocv(
        self, grids: np.ndarray, level: int
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """计算 pi 电子密度（MOCV）

        Args:
            grids: 三维坐标点列表 [[x, y, z], ...]
            level: 计算级别（0=仅密度，1=密度+梯度，2=密度+梯度+Hessian）

        Returns:
            (密度值, 梯度, Hessian) 三元组
        """
        return self.caler.pi_rho_mocv(grids, level)

    def ato_rho(
        self, grids: np.ndarray, level: int
    ) -> list[tuple[np.ndarray, np.ndarray, np.ndarray]]:
        """计算所有原子轨道在给定格点上的电子密度

        Args:
            grids: 三维坐标点列表 [[x, y, z], ...]
            level: 计算级别（0=仅密度，1=密度+梯度，2=密度+梯度+Hessian）

        Returns:
            每个原子轨道对应的 (密度值, 梯度, Hessian) 三元组列表
        """
        return self.caler.ato_rho(grids, level)

    def pro_rho(self, grids: np.ndarray) -> list[float]:
        """计算投影电子密度（仅含部分原子轨道贡献）

        Args:
            grids: 三维坐标点列表 [[x, y, z], ...]

        Returns:
            每个格点上的投影密度值
        """
        return self.caler.pro_rho(grids)

    def obt_rho(
        self, grids: np.ndarray, obt: int, level: int
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """计算指定分子轨道在给定格点上的电子密度

        Args:
            grids: 三维坐标点列表 [[x, y, z], ...]
            obt: 分子轨道索引（从 0 开始）
            level: 计算级别（0=仅密度，1=密度+梯度，2=密度+梯度+Hessian）

        Returns:
            (密度值, 梯度, Hessian) 三元组
        """
        return self.caler.obt_rho(grids, obt, level)

    def RDG(self, grids: np.ndarray) -> list[float]:
        """计算约化密度梯度（Reduced Density Gradient）

        用于描述弱相互作用区域，是 RDG / NCI 分析的核心量。

        Args:
            grids: 三维坐标点列表 [[x, y, z], ...]

        Returns:
            每个格点上的 RDG 值
        """
        return self.caler.RDG(grids)

    def IRI(self, grids: np.ndarray, alp: float | None = None) -> list[float]:
        """计算 IRI（Interaction Region Indicator）值

        IRI 用于可视化分子内和分子间的弱相互作用区域。

        Args:
            grids: 三维坐标点列表 [[x, y, z], ...]
            alp: 可选参数，指数缩放因子；若为 None 则使用默认值

        Returns:
            每个格点上的 IRI 值
        """
        return self.caler.IRI(grids, alp)

    def laplace(self, grids: np.ndarray) -> list[float]:
        """计算电子密度的拉普拉斯算子

        Args:
            grids: 三维坐标点列表 [[x, y, z], ...]

        Returns:
            每个格点上的拉普拉斯值
        """
        return self.caler.laplace(grids)

    def signl2rho(self, grids: np.ndarray) -> list[float]:
        """计算 sign(lambda2)*rho 值

        该量用于 NCI（非共价相互作用）分析，符号由 Hessian 矩阵的
        第二特征值 lambda2 决定。

        Args:
            grids: 三维坐标点列表 [[x, y, z], ...]

        Returns:
            每个格点上的 sign(lambda2)*rho 值
        """
        return self.caler.signl2rho(grids)

    def del_rho(self, grids: np.ndarray) -> list[float]:
        """计算电子密度的梯度模长 ||∇ρ||

        Args:
            grids: 三维坐标点列表 [[x, y, z], ...]

        Returns:
            每个格点上的梯度模长
        """
        return self.caler.del_rho(grids)

    def fukui(self, grids: np.ndarray, moleN: Mole, moleP: Mole) -> list[list[float]]:
        """计算 Fukui 函数值（亲核、亲电、自由基攻击）

        利用 N-1 和 N+1 电子态的分子计算 Fukui 函数。

        Args:
            grids: 三维坐标点列表 [[x, y, z], ...]
            moleN: N-1 电子态的分子对象
            moleP: N+1 电子态的分子对象

        Returns:
            每个格点上的 [f+, f-, f0, ...] 值列表
        """
        return self.caler.fukui(grids, moleN.core, moleP.core)

    def lagkin(self, grids: np.ndarray) -> list[list[float]]:
        """计算 Lagrangian 动能密度

        Args:
            grids: 三维坐标点列表 [[x, y, z], ...]

        Returns:
            每个格点上的 [K, Kxx, Kyy, Kzz] 值列表
        """
        return self.caler.lagkin(grids)

    def hamkin(self, grids: np.ndarray) -> list[float]:
        """计算 Hamiltonian 动能密度

        Args:
            grids: 三维坐标点列表 [[x, y, z], ...]

        Returns:
            每个格点上的 Hamiltonian 动能密度值
        """
        return self.caler.hamkin(grids)

    def KED(self, grids: np.ndarray) -> list[float]:
        """计算动能密度（Kinetic Energy Density）

        Args:
            grids: 三维坐标点列表 [[x, y, z], ...]

        Returns:
            每个格点上的动能密度值
        """
        return self.caler.KED(grids)

    def ELF(self, grids: np.ndarray) -> list[float]:
        """计算电子局域化函数（Electron Localization Function）

        ELF 用于描述电子的局域化程度，取值范围通常为 [0, 1]。

        Args:
            grids: 三维坐标点列表 [[x, y, z], ...]

        Returns:
            每个格点上的 ELF 值
        """
        return self.caler.ELF(grids)

    def LOL(self, grids: np.ndarray) -> list[float]:
        """计算局域化轨道指示函数（Localized Orbital Locator）

        LOL 与 ELF 类似，用于表征电子的局域化特征。

        Args:
            grids: 三维坐标点列表 [[x, y, z], ...]

        Returns:
            每个格点上的 LOL 值
        """
        return self.caler.LOL(grids)

    def infoentro(self, grids: np.ndarray) -> list[float]:
        """计算信息熵密度（Information Entropy Density）

        Args:
            grids: 三维坐标点列表 [[x, y, z], ...]

        Returns:
            每个格点上的信息熵密度值
        """
        return self.caler.infoentro(grids)

    def mole_pi_pocv(self) -> Mole:
        """获取 π 电子（POCV 方法）对应的分子对象"""
        return Mole(self.caler.mole_pi_pocv())

    def mole_pi_mocv(self) -> Mole:
        """获取 π 电子（MOCV 方法）对应的分子对象"""
        return Mole(self.caler.mole_pi_mocv())

    def mole_spin(self) -> tuple[Mole, Mole]:
        """获取自旋密度对应的两个分子对象（alpha, beta）"""
        a, b = self.caler.mole_spin()
        return Mole(a), Mole(b)

    def mol_rho_spin(self, grid: list[float], level: int = 0):
        """计算自旋密度

        Args:
            grid: 格点坐标
            level: 计算级别

        Returns:
            (alpha 密度, beta 密度)，均为 `GridVals` 对象
        """
        from pywfn.props.gridprop import GridVals

        a, b = self.caler.mol_rho_spin(grid, level)
        return GridVals(a), GridVals(b)

    def activate(self, grid: list[float]) -> float:
        """计算指定格点上的活化能指标"""
        return self.caler.activate(grid)

    def curvature(self, grid: list[float]) -> float:
        """计算指定格点上的曲率"""
        return self.caler.curvature(grid)

    def nuc_esp(self, grid: list[float]) -> float:
        """计算指定格点上的核静电势"""
        return self.caler.nuc_esp(grid)

    def ele_esp(self, grid: list[float]) -> float:
        """计算指定格点上的电子静电势"""
        return self.caler.ele_esp(grid)

    def tot_esp(self, grid: list[float]) -> float:
        """计算指定格点上的总静电势"""
        return self.caler.tot_esp(grid)
