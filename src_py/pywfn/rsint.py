"""
电子积分：单电子积分（重叠、动能、势能）与双电子积分

对应 rswfn 的 `rsint` 模块，Rust 侧由 `_core.rsint` 提供。
"""

import numpy as np

from pywfn import _core
from pywfn.base import Mole, Stm


class RsInt:
    """单电子积分（重叠、动能、势能）计算器"""

    def __init__(self, mole: Mole) -> None:
        """单电子积分计算器

        Args:
            mole: 分子对象
        """
        self.mole = mole
        self.core = _core.rsint.RsInt(mole.core)  # type: ignore

    def smat(self, form: str) -> np.ndarray:
        """计算重叠矩阵

        Args:
            form: 积分格式，"sph"（球谐）或 "car"（笛卡尔）

        Returns:
            np.ndarray: 重叠矩阵
        """
        return self.core.smat(form)

    def tmat(self, form: str) -> np.ndarray:
        """计算动能矩阵

        Args:
            form: 积分格式，"sph" 或 "car"

        Returns:
            np.ndarray: 动能矩阵
        """
        return self.core.tmat(form)

    def vmat(self, form: str) -> np.ndarray:
        """计算势能矩阵

        Args:
            form: 积分格式，"sph" 或 "car"

        Returns:
            np.ndarray: 势能（核吸引）矩阵
        """
        return self.core.vmat(form)

    def mat1e(self, cals: list[int], form: str) -> list[np.ndarray]:
        """计算单电子积分矩阵

        Args:
            cals: 指定需要计算的矩阵，[calc_s, calc_t, calc_v]，1 表示计算，0 表示跳过
            form: 积分格式，"sph" 或 "car"

        Returns:
            list[np.ndarray]: 三个二维数组 [S, T, V]
        """
        return self.core.mat1e(cals, form)

    def mat1e_sph(self, car_mat: np.ndarray) -> np.ndarray:
        """把笛卡尔基下的单电子积分矩阵变换到球谐基

        Args:
            car_mat: 笛卡尔基下的积分矩阵

        Returns:
            np.ndarray: 球谐基下的积分矩阵
        """
        return self.core.mat1e_sph(car_mat)

    def mat1e_raw(self, car_mat: np.ndarray) -> np.ndarray:
        """把积分矩阵变换回未归一化的原始基

        Args:
            car_mat: 笛卡尔基下的积分矩阵

        Returns:
            np.ndarray: 原始基下的积分矩阵
        """
        return self.core.mat1e_raw(car_mat)

    def mat1e_lcs(self, stms: dict[int, Stm], cals: list[int]) -> list[np.ndarray]:
        """计算局部坐标系下的单电子积分矩阵

        Args:
            stms: {原子索引: 局域坐标系}
            cals: 指定需要计算的矩阵，[calc_s, calc_t, calc_v]

        Returns:
            list[np.ndarray]: 三个二维数组 [S, T, V]
        """
        stms = {atm: stm.core for atm, stm in stms.items()}
        return self.core.mat1e_lcs(stms, cals)


class Mat2e:
    """双电子积分张量

    底层按四个基函数下标 `(i, j, k, l)` 压缩存储积分值。
    """

    def __init__(self, nato: int) -> None:
        """创建一个元素全为零的双电子积分张量

        Args:
            nato: 原子轨道（基函数）数量
        """
        self.core = _core.rsint.Mat2e(nato)  # type: ignore

    @staticmethod
    def from_core(core: "_core.rsint.Mat2e") -> "Mat2e":
        """把 `_core` 的 Mat2e 包装成 Python 对象

        Args:
            core: `_core.rsint.Mat2e` 对象

        Returns:
            Mat2e: 包装后的积分张量
        """
        mat = Mat2e.__new__(Mat2e)
        mat.core = core
        return mat

    @property
    def nato(self) -> int:
        """原子轨道（基函数）数量"""
        return self.core.nato

    @nato.setter
    def nato(self, nato: int) -> None:
        self.core.nato = nato

    @property
    def vals(self) -> list[float]:
        """一维压缩存储的积分值"""
        return self.core.vals

    @vals.setter
    def vals(self, vals: list[float]) -> None:
        self.core.vals = vals

    def get(self, i: int, j: int, k: int, l: int) -> float:
        """按四个基函数下标读取积分值

        Args:
            i: 第一个基函数下标
            j: 第二个基函数下标
            k: 第三个基函数下标
            l: 第四个基函数下标

        Returns:
            float: 积分值 (ij|kl)
        """
        return self.core.get(i, j, k, l)

    def set(self, i: int, j: int, k: int, l: int, val: float) -> None:
        """按四个基函数下标写入积分值

        Args:
            i: 第一个基函数下标
            j: 第二个基函数下标
            k: 第三个基函数下标
            l: 第四个基函数下标
            val: 积分值 (ij|kl)
        """
        self.core.set(i, j, k, l, val)

    def get_by_pair(self, p: int, q: int) -> float:
        """按壳层对下标读取积分值

        Args:
            p: 第一个壳层对下标
            q: 第二个壳层对下标

        Returns:
            float: 该壳层对对应的积分值
        """
        return self.core.get_by_pair(p, q)


class Engine:
    """双电子积分引擎

    负责双电子积分的预筛与分块。创建时会构建一次引擎并记录构建耗时、
    分块数、内存占用等信息。
    """

    def __init__(self, mole: Mole, form: str, screen: float, max_mem: int) -> None:
        """创建双电子积分引擎

        Args:
            mole: 分子对象
            form: 基函数形式，"sph" 或 "car"
            screen: 积分预筛阈值
            max_mem: 允许使用的最大内存（字节）
        """
        self.core = _core.rsint.Engine(mole.core, form, screen, max_mem)  # type: ignore

    @property
    def mole(self) -> Mole:
        """引擎所属的分子对象"""
        return Mole(self.core.mole)

    @mole.setter
    def mole(self, mole: Mole) -> None:
        self.core.mole = mole.core

    @property
    def form(self) -> str:
        """基函数形式，"sph" 或 "car" """
        return self.core.form

    @form.setter
    def form(self, form: str) -> None:
        self.core.form = form

    @property
    def screen(self) -> float:
        """积分预筛阈值"""
        return self.core.screen

    @screen.setter
    def screen(self, screen: float) -> None:
        self.core.screen = screen

    @property
    def max_mem(self) -> int:
        """允许使用的最大内存（字节）"""
        return self.core.max_mem

    @max_mem.setter
    def max_mem(self, max_mem: int) -> None:
        self.core.max_mem = max_mem

    @property
    def time_build(self) -> float:
        """引擎构建耗时（秒）"""
        return self.core.time_build

    @property
    def nblock(self) -> int:
        """积分分块数量"""
        return self.core.nblock

    @property
    def mem_bytes(self) -> int:
        """引擎预估的内存占用（字节）"""
        return self.core.mem_bytes

    def debug_check_blocks(self) -> tuple[int, int, float]:
        """检查积分分块信息

        Returns:
            tuple[int, int, float]: (分块数量, 原子轨道数量, 平均每个分块的积分数量)
        """
        return self.core.debug_check_blocks()

    def nato(self) -> int:
        """原子轨道（基函数）数量

        Returns:
            int: 基函数个数
        """
        return self.core.nato()

    def is_direct(self) -> bool:
        """是否为直接（direct）积分算法

        Returns:
            bool: 直接算法返回 True，否则返回 False
        """
        return self.core.is_direct()

    def mat2e_local_pi(
        self,
        tmats: list[np.ndarray],
        pi_comp: list[list[tuple[int, int]]],
        npi: int,
    ) -> Mat2e:
        """在局域坐标系下构建双电子积分张量

        Args:
            tmats: 各原子的局域坐标变换矩阵
            pi_comp: 每个原子的壳层对索引列表
            npi: 壳层对数量

        Returns:
            Mat2e: 构建好的双电子积分张量
        """
        return Mat2e.from_core(self.core.mat2e_local_pi(tmats, pi_comp, npi))

    def jk(self, dmat: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """计算库仑（J）与交换（K）矩阵

        Args:
            dmat: 密度矩阵

        Returns:
            tuple[np.ndarray, np.ndarray]: (J 矩阵, K 矩阵)
        """
        return self.core.jk(dmat)

    def jk_spin(
        self, d_tot: np.ndarray, d_a: np.ndarray, d_b: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """计算开壳层的库仑（J）与自旋分辨交换（K）矩阵

        Args:
            d_tot: 总密度矩阵
            d_a: alpha 自旋密度矩阵
            d_b: beta 自旋密度矩阵

        Returns:
            tuple[np.ndarray, np.ndarray, np.ndarray]: (J, K_alpha, K_beta)
        """
        return self.core.jk_spin(d_tot, d_a, d_b)


class EspCalculator:
    """静电势（ESP）计算器"""

    def __init__(self, mole: Mole, dmat: np.ndarray) -> None:
        """创建静电势计算器

        Args:
            mole: 分子对象
            dmat: 密度矩阵
        """
        self.core = _core.rsint.EspCalculator(mole.core, dmat)  # type: ignore

    @property
    def mole(self) -> Mole:
        """计算所用的分子对象"""
        return Mole(self.core.mole)

    @property
    def dmat(self) -> np.ndarray:
        """计算所用的密度矩阵"""
        return np.array(self.core.dmat)

    def val(self, rc: list[float]) -> float:
        """计算单个点的静电势

        Args:
            rc: 点的笛卡尔坐标 [x, y, z]

        Returns:
            float: 该点的静电势值
        """
        return self.core.val(rc)

    def vals(self, rcs: list[list[float]]) -> list[float]:
        """批量计算多个点的静电势

        Args:
            rcs: 点坐标列表，每个元素为 [x, y, z]

        Returns:
            list[float]: 各点的静电势值
        """
        return self.core.vals(rcs)


def get_2eid(i: int, j: int, k: int, l: int) -> int:
    """计算双电子积分在一维压缩存储中的下标

    Args:
        i: 第一个基函数下标
        j: 第二个基函数下标
        k: 第三个基函数下标
        l: 第四个基函数下标

    Returns:
        int: 对应的一维下标
    """
    return _core.rsint.get_2eid(i, j, k, l)  # type: ignore


def jk_mat(dmat: np.ndarray, emat: Mat2e) -> tuple[np.ndarray, np.ndarray]:
    """由积分张量与密度矩阵计算库仑（J）与交换（K）矩阵

    Args:
        dmat: 密度矩阵
        emat: 双电子积分张量

    Returns:
        tuple[np.ndarray, np.ndarray]: (J 矩阵, K 矩阵)
    """
    return _core.rsint.jk_mat(dmat, emat.core)  # type: ignore
