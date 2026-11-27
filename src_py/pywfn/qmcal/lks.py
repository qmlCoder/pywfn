"""
局部 π 电子 Kohn-Sham（LKS）计算

对应 rswfn 的 `qmcal::lks` 模块，只对分子的 π 电子做 Kohn-Sham 自洽场计算，
σ 骨架以冻结的方式提供势场。
"""

from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray

from pywfn import _core
from pywfn.base.mole import Mole
from pywfn.base.atoms import Stm


@dataclass
class LKSResult:
    """局部 π 电子 Kohn-Sham（LKS）计算结果"""

    ele_eng: float
    xc_eng: float
    nuc_eng: float
    converged: bool
    niter: int
    delta_eng: float
    npi_elec: int
    """π 电子数"""
    npi: int
    """π 轨道数"""
    nato: int
    """原子轨道数"""
    ngrid: int
    obt_engs: list[float]
    pi_atos: list[int]
    cmat: np.ndarray
    """π 轨道系数矩阵（nato × npi）"""
    _core: "_core.qmcal.lks.LKSResult | None" = field(
        default=None, repr=False, compare=False
    )
    """底层 `_core` 结果对象，供展开系数、生成轨道表格等方法使用"""

    @property
    def total_eng(self) -> float:
        """分子总能量（电子能量 + 核排斥能）"""
        return self.ele_eng + self.nuc_eng

    @property
    def nocc(self) -> int:
        """占据轨道数"""
        return (self.npi_elec + 1) // 2

    def _res(self) -> "_core.qmcal.lks.LKSResult":
        """获取底层计算结果对象

        Returns:
            _core.qmcal.lks.LKSResult: 底层结果对象

        Raises:
            RuntimeError: 该结果对象不是由 `LKS.run` 生成时
        """
        if self._core is None:
            raise RuntimeError("该结果对象未绑定底层计算结果，请用 LKS.run 生成")
        return self._core

    @property
    def func(self) -> str:
        """泛函名称"""
        return self._res().func

    @property
    def stms(self) -> dict[int, Stm]:
        """各 π 原子轨道对应的局域坐标系"""
        return {atm: Stm(stm) for atm, stm in self._res().stms.items()}

    @property
    def time_int1e(self) -> float:
        """单电子积分耗时（秒）"""
        return self._res().time_int1e

    @property
    def time_int2e(self) -> float:
        """双电子积分耗时（秒）"""
        return self._res().time_int2e

    @property
    def time_scf(self) -> float:
        """自洽场迭代总耗时（秒）"""
        return self._res().time_scf

    def expand_coefs(self, mole: Mole) -> NDArray[np.float64]:
        """把局域 π 轨道系数展开为全基组的系数矩阵

        Args:
            mole: 分子对象

        Returns:
            NDArray[np.float64]: 展开后的系数矩阵
        """
        return self._res().expand_coefs(mole.core)

    def local_labels(self, mole: Mole) -> list[str]:
        """获取各局域 π 轨道的标签

        Args:
            mole: 分子对象

        Returns:
            list[str]: 轨道标签列表
        """
        return self._res().local_labels(mole.core)

    def mo_table_local(self, mole: Mole) -> str:
        """获取局域 π 轨道的表格化文本

        Args:
            mole: 分子对象

        Returns:
            str: 轨道信息表格
        """
        return self._res().mo_table_local(mole.core)

    def to_mole(self, mole: Mole) -> Mole:
        """把局域 π 计算结果还原成分子对象（写入轨道系数）

        Args:
            mole: 分子对象

        Returns:
            Mole: 带有本次计算轨道系数的新分子对象
        """
        return Mole(self._res().to_mole(mole.core))


class LKS:
    """局部 π 电子 Kohn-Sham（LKS）计算器

    Args:
        mole: 分子对象，需带有基组信息
        func: 泛函名，如 "b3lyp"
    """

    def __init__(self, mole: Mole, func: str) -> None:
        self.mole = mole
        self.core = _core.qmcal.lks.LKS(mole.core, func)  # type: ignore

    @classmethod
    def with_grid(
        cls,
        mole: Mole,
        radial: int,
        angular: int,
        func: str,
        npi_elec: int | None = None,
    ) -> "LKS":
        """用指定格点与 π 电子数创建计算器

        Args:
            mole: 分子对象
            radial: 径向格点数
            angular: 球面格点数
            func: 泛函名
            npi_elec: π 电子数；为 None 时由分子自动推断

        Returns:
            LKS: 局部 π 电子 KS 计算器
        """
        obj = cls.__new__(cls)
        obj.mole = mole
        obj.core = _core.qmcal.lks.LKS.with_grid(  # type: ignore
            mole.core, radial, angular, func, npi_elec
        )
        return obj

    @classmethod
    def from_name(cls, mole: Mole, name: str) -> "LKS | None":
        """按泛函名创建计算器（默认格点 100×74）

        Args:
            mole: 分子对象
            name: 泛函名，如 "b3lyp"

        Returns:
            LKS | None: 泛函名无效时返回 None
        """
        core = _core.qmcal.lks.LKS.from_name(mole.core, name)  # type: ignore
        if core is None:
            return None
        obj = cls.__new__(cls)
        obj.mole = mole
        obj.core = core
        return obj

    def run(self) -> LKSResult:
        """执行自洽场计算

        Returns:
            LKSResult: 计算结果
        """
        r = self.core.run()
        return LKSResult(
            ele_eng=r.ele_eng,
            xc_eng=r.xc_eng,
            nuc_eng=r.nuc_eng,
            converged=r.converged,
            niter=r.niter,
            delta_eng=r.delta_eng,
            npi_elec=r.npi_elec,
            npi=r.npi,
            nato=r.nato,
            ngrid=r.ngrid,
            obt_engs=r.obt_engs,
            pi_atos=r.pi_atos,
            cmat=r.cmat(),
            _core=r,
        )

    def nuc_eng(self) -> float:
        """计算核-核排斥能

        Returns:
            float: 核-核排斥能（Hartree）
        """
        return self.core.nuc_eng()

    @property
    def func(self) -> str:
        """泛函名称"""
        return self.core.func

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
    def npi_elec(self) -> int | None:
        """π 电子数；None 表示由分子自动推断"""
        return self.core.npi_elec

    @npi_elec.setter
    def npi_elec(self, value: int | None) -> None:
        self.core.npi_elec = value

    @property
    def stms(self) -> dict[int, Stm]:
        """各 π 原子轨道对应的局域坐标系"""
        return {atm: Stm(stm) for atm, stm in self.core.stms.items()}

    @property
    def pi_atos(self) -> list[int]:
        """π 原子轨道索引"""
        return self.core.pi_atos

    @property
    def grids(self) -> NDArray[np.float64]:
        """数值积分格点坐标，形状为 (ngrid, 3)"""
        return self.core.grids

    @property
    def weits(self) -> NDArray[np.float64]:
        """数值积分格点权重"""
        return self.core.weits

    @property
    def npi(self) -> int:
        """π 轨道数"""
        return self.core.npi

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
    def ts(self) -> list[NDArray[np.float64]]:
        """π 基函数之间的变换矩阵列表"""
        return self.core.ts

    @property
    def pi_comp(self) -> list[list[tuple[int, int]]]:
        """各 π 片段包含的原子轨道索引对"""
        return self.core.pi_comp

    @property
    def chi0(self) -> NDArray[np.float64]:
        """原子轨道在格点上的取值"""
        return self.core.chi0

    @property
    def chi1(self) -> NDArray[np.float64]:
        """原子轨道在格点上的梯度分量"""
        return self.core.chi1
