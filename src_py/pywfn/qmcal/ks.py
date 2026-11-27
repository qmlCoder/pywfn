"""
Kohn-Sham（DFT）计算

对应 rswfn 的 `qmcal::ks` 模块，提供闭壳层（RKS）与自旋极化（UKS）
Kohn-Sham 自洽场计算、数值积分格点相关数据结构，以及由矩阵计算电子能量的工具函数。
"""

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from pywfn import _core
from pywfn.base.mole import Mole

B3LYP_EXACT_X: float = _core.qmcal.ks.B3LYP_EXACT_X  # type: ignore
"""B3LYP 泛函中精确交换（Hartree-Fock 交换）的混合系数"""


@dataclass
class RKSResult:
    """闭壳层 Kohn-Sham（DFT）计算结果"""

    ele_eng: float
    """电子能量（含交换相关能）"""
    xc_eng: float
    """交换相关能"""
    nuc_eng: float
    """核-核排斥能"""
    converged: bool
    """是否自洽收敛"""
    niter: int
    """实际迭代次数"""
    delta_eng: float
    """最后一次迭代的能量变化"""
    form: str
    """实际使用的积分形式（"sph" 或 "car"）"""
    ngrid: int
    """数值积分格点数量"""

    @property
    def total_eng(self) -> float:
        """分子总能量（电子能量 + 核排斥能）"""
        return self.ele_eng + self.nuc_eng


@dataclass
class UKSResult:
    """自旋极化（开壳层）Kohn-Sham（DFT）计算结果"""

    ele_eng: float
    xc_eng: float
    nuc_eng: float
    converged: bool
    niter: int
    delta_eng: float
    form: str
    ngrid: int
    nela: int
    """α 电子数"""
    nelb: int
    """β 电子数"""
    spin_sq: float
    """⟨S²⟩ 期望值"""

    @property
    def total_eng(self) -> float:
        """分子总能量（电子能量 + 核排斥能）"""
        return self.ele_eng + self.nuc_eng

    @property
    def multip(self) -> int:
        """自旋多重度 2S + 1"""
        return self.nela - self.nelb + 1


class GridBlock:
    """数值积分格点分块信息

    描述数值积分格点中属于同一个原子的一段连续格点，用于按原子分批
    计算交换相关泛函。

    Args:
        p0: 分块起始格点索引
        p1: 分块结束格点索引
        atm: 该分块所属的原子索引
        rmin: 分块内离原子核的最近距离
        rmax: 分块内离原子核的最远距离
    """

    def __init__(self, p0: int, p1: int, atm: int, rmin: float, rmax: float) -> None:
        self.core = _core.qmcal.ks.GridBlock(p0, p1, atm, rmin, rmax)  # type: ignore

    @property
    def p0(self) -> int:
        """分块起始格点索引"""
        return self.core.p0

    @property
    def p1(self) -> int:
        """分块结束格点索引"""
        return self.core.p1

    @property
    def atm(self) -> int:
        """该分块所属的原子索引"""
        return self.core.atm

    @property
    def rmin(self) -> float:
        """分块内离原子核的最近距离"""
        return self.core.rmin

    @property
    def rmax(self) -> float:
        """分块内离原子核的最远距离"""
        return self.core.rmax

    def __repr__(self) -> str:
        return self.core.__repr__()


class SpinRho:
    """自旋极化电子密度

    保存 α/β 电子密度、它们对坐标的梯度，以及梯度的各种内积
    （σ_aa、σ_bb、σ_ab），供 GGA 类泛函使用。

    Args:
        core: 底层 `_core` 的 SpinRho 对象
    """

    def __init__(self, core: "_core.qmcal.ks.SpinRho") -> None:
        self.core = core

    @property
    def rho_a(self) -> list[float]:
        """α 电子密度"""
        return self.core.rho_a

    @property
    def rho_b(self) -> list[float]:
        """β 电子密度"""
        return self.core.rho_b

    @property
    def grad_a(self) -> NDArray[np.float64]:
        """α 密度梯度，形状为 (ngrid, 3)"""
        return self.core.grad_a

    @property
    def grad_b(self) -> NDArray[np.float64]:
        """β 密度梯度，形状为 (ngrid, 3)"""
        return self.core.grad_b

    @property
    def saa(self) -> list[float]:
        """α 密度梯度内积 `∇ρ_α · ∇ρ_α`"""
        return self.core.saa

    @property
    def sbb(self) -> list[float]:
        """β 密度梯度内积 `∇ρ_β · ∇ρ_β`"""
        return self.core.sbb

    @property
    def sab(self) -> list[float]:
        """α/β 密度梯度交叉内积 `∇ρ_α · ∇ρ_β`"""
        return self.core.sab

    def __repr__(self) -> str:
        return self.core.__repr__()


class RKS:
    """闭壳层 Kohn-Sham（DFT）计算器

    Args:
        mole: 分子对象，需带有基组信息
        func: 泛函名，如 "b3lyp"、"pbe0"，可由 `pywfn.rsxcs.supported_names` 查询
    """

    def __init__(self, mole: Mole, func: str) -> None:
        self.mole = mole
        self.core = _core.qmcal.ks.RKS(mole.core, func)  # type: ignore

    @classmethod
    def with_grid(cls, mole: Mole, radial: int, angular: int, func: str) -> "RKS":
        """用指定的径向/球面格点数创建计算器

        Args:
            mole: 分子对象
            radial: 径向格点数
            angular: 球面格点数
            func: 泛函名

        Returns:
            RKS: 闭壳层 KS 计算器
        """
        obj = cls.__new__(cls)
        obj.mole = mole
        obj.core = _core.qmcal.ks.RKS.with_grid(mole.core, radial, angular, func)  # type: ignore
        return obj

    @classmethod
    def with_grid_prune(
        cls, mole: Mole, radial: int, angular: int, func: str, prune: bool
    ) -> "RKS":
        """用指定的格点数创建计算器，并指定是否裁剪格点

        Args:
            mole: 分子对象
            radial: 径向格点数
            angular: 球面格点数
            func: 泛函名
            prune: 是否对格点做裁剪

        Returns:
            RKS: 闭壳层 KS 计算器
        """
        obj = cls.__new__(cls)
        obj.mole = mole
        obj.core = _core.qmcal.ks.RKS.with_grid_prune(  # type: ignore
            mole.core, radial, angular, func, prune
        )
        return obj

    @classmethod
    def from_name(cls, mole: Mole, name: str) -> "RKS | None":
        """按泛函名创建计算器（默认格点）

        Args:
            mole: 分子对象
            name: 泛函名，如 "b3lyp"

        Returns:
            RKS | None: 泛函名无效时返回 None
        """
        core = _core.qmcal.ks.RKS.from_name(mole.core, name)  # type: ignore
        if core is None:
            return None
        obj = cls.__new__(cls)
        obj.mole = mole
        obj.core = core
        return obj

    def run(self) -> RKSResult:
        """执行自洽场计算

        Returns:
            RKSResult: 计算结果
        """
        r = self.core.run()
        return RKSResult(
            ele_eng=r.ele_eng,
            xc_eng=r.xc_eng,
            nuc_eng=r.nuc_eng,
            converged=r.converged,
            niter=r.niter,
            delta_eng=r.delta_eng,
            form=r.form,
            ngrid=r.ngrid,
        )

    def nuc_eng(self) -> float:
        """计算核-核排斥能

        Returns:
            float: 核-核排斥能（Hartree）
        """
        return self.core.nuc_eng()

    def ao_stats(self) -> tuple[int, float, int]:
        """原子轨道统计信息

        Returns:
            tuple[int, float, int]: `(屏幕后有效原子轨道数, 最大系数绝对值, 原始原子轨道数)`
        """
        return self.core.ao_stats()

    def calc_rho(
        self, dmat: NDArray[np.float64]
    ) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
        """在数值积分格点上计算电子密度及其梯度

        Args:
            dmat: 密度矩阵，形状为 (nato, nato)

        Returns:
            tuple: `(格点密度, 格点梯度, 格点密度梯度平方)`，梯度形状为 (ngrid, 3)
        """
        return self.core.calc_rho(dmat)

    def xc_eng(self, dmat: NDArray[np.float64]) -> float:
        """计算交换相关能

        Args:
            dmat: 密度矩阵，形状为 (nato, nato)

        Returns:
            float: 交换相关能（Hartree）
        """
        return self.core.xc_eng(dmat)

    def grid_ele(self, dmat: NDArray[np.float64]) -> float:
        """计算电子在数值积分格点上的静电相互作用能

        Args:
            dmat: 密度矩阵，形状为 (nato, nato)

        Returns:
            float: 格点贡献的电子能量（Hartree）
        """
        return self.core.grid_ele(dmat)

    @property
    def func(self) -> str:
        """泛函名称"""
        return self.core.func

    @property
    def int_form(self) -> str:
        """积分形式（"sph" 或 "car"）"""
        return self.core.int_form

    @int_form.setter
    def int_form(self, value: str) -> None:
        self.core.int_form = value

    @property
    def radial(self) -> int:
        """径向格点数"""
        return self.core.radial

    @property
    def angular(self) -> int:
        """球面格点数"""
        return self.core.angular

    @property
    def max_iter(self) -> int:
        """最大 SCF 迭代次数"""
        return self.core.max_iter

    @max_iter.setter
    def max_iter(self, value: int) -> None:
        self.core.max_iter = value

    @property
    def tol(self) -> float:
        """SCF 收敛阈值"""
        return self.core.tol

    @tol.setter
    def tol(self, value: float) -> None:
        self.core.tol = value

    @property
    def grids(self) -> NDArray[np.float64]:
        """数值积分格点坐标，形状为 (ngrid, 3)"""
        return self.core.grids

    @property
    def weits(self) -> NDArray[np.float64]:
        """数值积分格点权重"""
        return self.core.weits

    @property
    def ngrid(self) -> int:
        """数值积分格点数量"""
        return self.core.ngrid

    @property
    def nato(self) -> int:
        """原子轨道数量"""
        return self.core.nato

    @property
    def chi0(self) -> NDArray[np.float64]:
        """原子轨道在格点上的取值"""
        return self.core.chi0

    @property
    def chi1(self) -> NDArray[np.float64]:
        """原子轨道在格点上的梯度分量"""
        return self.core.chi1

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

    @property
    def ao_cut(self) -> float:
        """原子轨道在格点上的截断阈值"""
        return self.core.ao_cut

    @ao_cut.setter
    def ao_cut(self, value: float) -> None:
        self.core.ao_cut = value


class UKS:
    """自旋极化（开壳层）Kohn-Sham（DFT）计算器

    Args:
        mole: 分子对象，需带有基组信息
        func: 泛函名，如 "b3lyp"、"pbe0"，可由 `pywfn.rsxcs.supported_names` 查询
    """

    def __init__(self, mole: Mole, func: str) -> None:
        self.mole = mole
        self.core = _core.qmcal.ks.UKS(mole.core, func)  # type: ignore

    @classmethod
    def with_grid(cls, mole: Mole, radial: int, angular: int, func: str) -> "UKS":
        """用指定的径向/球面格点数创建计算器

        Args:
            mole: 分子对象
            radial: 径向格点数
            angular: 球面格点数
            func: 泛函名

        Returns:
            UKS: 自旋极化 KS 计算器
        """
        obj = cls.__new__(cls)
        obj.mole = mole
        obj.core = _core.qmcal.ks.UKS.with_grid(mole.core, radial, angular, func)  # type: ignore
        return obj

    @classmethod
    def with_grid_prune(
        cls, mole: Mole, radial: int, angular: int, func: str, prune: bool
    ) -> "UKS":
        """用指定的格点数创建计算器，并指定是否裁剪格点

        Args:
            mole: 分子对象
            radial: 径向格点数
            angular: 球面格点数
            func: 泛函名
            prune: 是否对格点做裁剪

        Returns:
            UKS: 自旋极化 KS 计算器
        """
        obj = cls.__new__(cls)
        obj.mole = mole
        obj.core = _core.qmcal.ks.UKS.with_grid_prune(  # type: ignore
            mole.core, radial, angular, func, prune
        )
        return obj

    @classmethod
    def from_name(cls, mole: Mole, name: str) -> "UKS | None":
        """按泛函名创建计算器（默认格点）

        Args:
            mole: 分子对象
            name: 泛函名，如 "b3lyp"

        Returns:
            UKS | None: 泛函名无效时返回 None
        """
        core = _core.qmcal.ks.UKS.from_name(mole.core, name)  # type: ignore
        if core is None:
            return None
        obj = cls.__new__(cls)
        obj.mole = mole
        obj.core = core
        return obj

    def run(self) -> UKSResult:
        """执行自洽场计算

        Returns:
            UKSResult: 计算结果
        """
        r = self.core.run()
        return UKSResult(
            ele_eng=r.ele_eng,
            xc_eng=r.xc_eng,
            nuc_eng=r.nuc_eng,
            converged=r.converged,
            niter=r.niter,
            delta_eng=r.delta_eng,
            form=r.form,
            ngrid=r.ngrid,
            nela=r.nela,
            nelb=r.nelb,
            spin_sq=r.spin_sq,
        )

    def nuc_eng(self) -> float:
        """计算核-核排斥能

        Returns:
            float: 核-核排斥能（Hartree）
        """
        return self.core.nuc_eng()

    def calc_rho_spin(
        self, dmat_a: NDArray[np.float64], dmat_b: NDArray[np.float64]
    ) -> SpinRho:
        """在数值积分格点上计算 α/β 自旋密度及其梯度

        Args:
            dmat_a: α 密度矩阵，形状为 (nato, nato)
            dmat_b: β 密度矩阵，形状为 (nato, nato)

        Returns:
            SpinRho: 自旋密度对象（含密度、梯度及各密度梯度内积）
        """
        return SpinRho(self.core.calc_rho_spin(dmat_a, dmat_b))

    @property
    def base(self) -> RKS:
        """底层闭壳层 KS 计算器"""
        obj = RKS.__new__(RKS)
        obj.mole = self.mole
        obj.core = self.core.base
        return obj

    @property
    def nela(self) -> int:
        """α 电子数"""
        return self.core.nela

    @property
    def nelb(self) -> int:
        """β 电子数"""
        return self.core.nelb


def rks_ele_eng(
    dmat: NDArray[np.float64],
    hmat: NDArray[np.float64],
    jmat: NDArray[np.float64],
    kmat: NDArray[np.float64],
    exact_x: float,
    xc_eng: float,
) -> float:
    """计算闭壳层 Kohn-Sham 电子能量

    Args:
        dmat: 密度矩阵，形状为 (nato, nato)
        hmat: 核心哈密顿矩阵，形状为 (nato, nato)
        jmat: 库仑矩阵，形状为 (nato, nato)
        kmat: 交换矩阵，形状为 (nato, nato)
        exact_x: 精确交换的混合系数
        xc_eng: 交换相关能（Hartree）

    Returns:
        float: 电子能量（Hartree）
    """
    return _core.qmcal.ks.rks_ele_eng(  # type: ignore
        dmat, hmat, jmat, kmat, exact_x, xc_eng
    )


def uks_ele_eng(
    dmat: NDArray[np.float64],
    hmat: NDArray[np.float64],
    jmat: NDArray[np.float64],
    dmat_a: NDArray[np.float64],
    dmat_b: NDArray[np.float64],
    kmat_a: NDArray[np.float64],
    kmat_b: NDArray[np.float64],
    exact_x: float,
    xc_eng: float,
) -> float:
    """计算自旋极化 Kohn-Sham 电子能量

    Args:
        dmat: 总密度矩阵，形状为 (nato, nato)
        hmat: 核心哈密顿矩阵，形状为 (nato, nato)
        jmat: 库仑矩阵，形状为 (nato, nato)
        dmat_a: α 密度矩阵，形状为 (nato, nato)
        dmat_b: β 密度矩阵，形状为 (nato, nato)
        kmat_a: α 交换矩阵，形状为 (nato, nato)
        kmat_b: β 交换矩阵，形状为 (nato, nato)
        exact_x: 精确交换的混合系数
        xc_eng: 交换相关能（Hartree）

    Returns:
        float: 电子能量（Hartree）
    """
    return _core.qmcal.ks.uks_ele_eng(  # type: ignore
        dmat, hmat, jmat, dmat_a, dmat_b, kmat_a, kmat_b, exact_x, xc_eng
    )
