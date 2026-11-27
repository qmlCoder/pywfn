"""
Hartree-Fock 计算

对应 rswfn 的 `qmcal::hf` 模块，提供闭壳层（RHF）与无限制（UHF）Hartree-Fock
自洽场计算，以及由密度矩阵与 Fock 矩阵计算电子能量的工具函数。
"""

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from pywfn import _core
from pywfn.base import Basis, Coefs, Mole


@dataclass
class RHFResult:
    """闭壳层（限制性）Hartree-Fock（RHF）计算结果"""

    ele_eng: float
    """电子能量（Hartree）"""
    nuc_eng: float
    """核-核排斥能（Hartree）"""
    converged: bool
    """是否自洽收敛"""
    niter: int
    """实际迭代次数"""
    delta_eng: float
    """最后一次迭代的能量变化（Hartree）"""
    form: str
    """实际使用的积分形式（"sph" 或 "car"）"""
    time_int1e: float
    """单电子积分耗时（秒）"""
    time_int2e: float
    """双电子积分耗时（秒）"""
    time_scf: float
    """自洽场迭代总耗时（秒）"""
    time_jk: float
    """构建 J/K 矩阵耗时（秒）"""
    coefs: Coefs
    """分子轨道系数"""
    basis: Basis
    """基组信息"""

    @property
    def total_eng(self) -> float:
        """分子总能量（电子能量 + 核排斥能）"""
        return self.ele_eng + self.nuc_eng


class RHF:
    """RHF（闭壳层 Hartree-Fock）计算器

    求解闭壳层 Hartree-Fock 方程，双电子积分由 rswfn 内部的壳层块引擎
    直接计算（直接 SCF 算法）。

    Args:
        mole: 分子对象，需带有基组信息
    """

    def __init__(self, mole: Mole) -> None:
        self.mole = mole
        self.core = _core.qmcal.hf.RHF(mole.core)  # type: ignore

    def run(self) -> RHFResult:
        """执行 RHF 自洽场计算

        迭代求解 Hartree-Fock 方程，获得收敛的分子轨道与能量。

        Returns:
            RHFResult: 计算结果，包含能量、迭代信息、轨道系数与基组
        """
        r = self.core.run()
        # Coefs 包装类的 get_cmat/get_dmat 统一从所绑定的分子读取系数，
        # 因此先用计算结果重建一个分子对象，再从中取 Coefs 与 Basis
        res_mole = Mole.build(self.mole.atoms, self.mole.basis, Coefs(r.coefs, self.mole))
        return RHFResult(
            ele_eng=r.ele_eng,
            nuc_eng=r.nuc_eng,
            converged=r.converged,
            niter=r.niter,
            delta_eng=r.delta_eng,
            form=r.form,
            time_int1e=r.time_int1e,
            time_int2e=r.time_int2e,
            time_scf=r.time_scf,
            time_jk=r.time_jk,
            coefs=res_mole.coefs,
            basis=Basis(r.basis, res_mole),
        )

    def nuc_eng(self) -> float:
        """计算核-核排斥能

        Returns:
            float: 核-核排斥能（Hartree）
        """
        return self.core.nuc_eng()

    @property
    def int_form(self) -> str:
        """积分形式（"sph" 或 "car"）"""
        return self.core.int_form

    @int_form.setter
    def int_form(self, value: str) -> None:
        self.core.int_form = value

    @property
    def max_iter(self) -> int:
        """最大 SCF 迭代次数"""
        return self.core.max_iter

    @max_iter.setter
    def max_iter(self, value: int) -> None:
        self.core.max_iter = value

    @property
    def tol(self) -> float:
        """SCF 收敛阈值（能量变化，Hartree）"""
        return self.core.tol

    @tol.setter
    def tol(self, value: float) -> None:
        self.core.tol = value

    @property
    def screen(self) -> float:
        """双电子积分 screening 阈值"""
        return self.core.screen

    @screen.setter
    def screen(self, value: float) -> None:
        self.core.screen = value

    @property
    def max_mem(self) -> int:
        """内存上限（字节）"""
        return self.core.max_mem

    @max_mem.setter
    def max_mem(self, value: int) -> None:
        self.core.max_mem = value


class UHF:
    """UHF（无限制 Hartree-Fock）计算器

    对 α/β 电子分别求解，适用于开壳层体系。

    Args:
        mole: 分子对象，需带有基组信息
    """

    def __init__(self, mole: Mole) -> None:
        self.mole = mole
        self.core = _core.qmcal.hf.UHF(mole.core)  # type: ignore

    def calc_udmat(
        self,
        nele: list[int] | tuple[int, int],
        cmat_a: NDArray[np.float64],
        cmat_b: NDArray[np.float64],
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """由 α/β 轨道系数构造密度矩阵

        Args:
            nele: α 与 β 电子数，形如 [nela, nelb]
            cmat_a: α 轨道系数矩阵，形状为 (nato, nato)
            cmat_b: β 轨道系数矩阵，形状为 (nato, nato)

        Returns:
            tuple[NDArray[np.float64], NDArray[np.float64]]: (α 密度矩阵, β 密度矩阵)
        """
        return self.core.calc_udmat(list(nele), cmat_a, cmat_b)

    def calc_uhf(
        self,
    ) -> tuple[
        list[float],
        tuple[NDArray[np.float64], NDArray[np.float64]],
        tuple[NDArray[np.float64], NDArray[np.float64]],
    ]:
        """执行无限制 Hartree-Fock 计算

        Returns:
            tuple: 三元组 `(α/β 电子能量, (α 轨道能量, β 轨道能量), (α 轨道系数, β 轨道系数))`
        """
        return self.core.calc_uhf()


def rhf_ele_eng(
    dmat: NDArray[np.float64],
    hmat: NDArray[np.float64],
    fmat: NDArray[np.float64],
) -> float:
    """计算闭壳层电子能量

    Args:
        dmat: 密度矩阵，形状为 (nato, nato)
        hmat: 核心哈密顿矩阵（动能 + 核吸引），形状为 (nato, nato)
        fmat: Fock 矩阵，形状为 (nato, nato)

    Returns:
        float: 电子能量（Hartree）
    """
    return _core.qmcal.hf.rhf_ele_eng(dmat, hmat, fmat)  # type: ignore


def uhf_ele_eng(
    dmat_a: NDArray[np.float64],
    dmat_b: NDArray[np.float64],
    tmat: NDArray[np.float64],
    vmat: NDArray[np.float64],
    fmat_a: NDArray[np.float64],
    fmat_b: NDArray[np.float64],
) -> list[float]:
    """计算开壳层（无限制）电子能量

    Args:
        dmat_a: α 密度矩阵，形状为 (nato, nato)
        dmat_b: β 密度矩阵，形状为 (nato, nato)
        tmat: 动能矩阵，形状为 (nato, nato)
        vmat: 核吸引矩阵，形状为 (nato, nato)
        fmat_a: α Fock 矩阵，形状为 (nato, nato)
        fmat_b: β Fock 矩阵，形状为 (nato, nato)

    Returns:
        list[float]: `[α 电子能量, β 电子能量]`（Hartree）
    """
    return _core.qmcal.hf.uhf_ele_eng(  # type: ignore
        dmat_a, dmat_b, tmat, vmat, fmat_a, fmat_b
    )
