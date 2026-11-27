from typing import Literal

import numpy as np

from pywfn import _core
from pywfn.base import Mole, Stm

Chrgs = Literal["mulliken", "lowdin", "space", "hirshfeld"]


class Calculator:
    def __init__(self, mole: "Mole", atms: list[int] | None = None):
        self.mole = mole
        self.caler = _core.props.atomprop.charge.Calculator(mole.core, atms)  # type: ignore # 核心计算器

    def spin(self, ctype: str = "mulliken") -> np.ndarray:
        """计算所有原子的自旋"""
        return self.caler.spin(ctype)

    def mulliken(self) -> np.ndarray:
        """计算Mulliken电荷"""
        return self.caler.mulliken()

    def lowdin(self) -> np.ndarray:
        """计算每个原子的Lowdin电荷"""
        return self.caler.lowdin()

    def hirshfeld(self) -> np.ndarray:
        """计算原子的Hirshfeld电荷"""
        return self.caler.hirshfeld()

    def pocv(
        self,
        dirs: dict[int, list[float]],
        keep_other_atm: bool,
        keep_other_sym: bool,
        ctype: str,
    ) -> np.ndarray:
        """计算投影分子轨道的电子数"""
        return self.caler.pocv(dirs, keep_other_atm, keep_other_sym, ctype)  # type: ignore

    def mocv(
        self, stms: dict[int, Stm], atos: list[int], ctype: str = "mulliken"
    ) -> np.ndarray:
        """计算MOCV电荷

        Args:
            stms: 原子局部坐标系
            atos: 原子轨道索引
            ctype: 电荷类型, 默认为 mulliken

        Returns:
            np.ndarray: MOCV电荷
        """
        stms = {atm: stm.core for atm, stm in stms.items()}
        return self.caler.mocv(stms, atos, ctype)

    def pi_pocv(
        self, ctype: str = "mulliken"
    ) -> tuple[dict[int, list[float]], np.ndarray]:
        """
        计算π电子数
        每个原子的方向为法向量的方向电子数即为π电子数
        """
        dirs, eles = self.caler.pi_pocv(ctype)
        return dirs, np.array(eles)

    def pi_mocv(
        self, ctype: str = "mulliken"
    ) -> tuple[dict[int, Stm], list[int], np.ndarray]:
        """使用轨道分解方法计算pi电子分布, 可以包含D轨道"""
        stms, atos, eles = self.caler.pi_mocv(ctype)
        stms = {atm: Stm(stm) for atm, stm in stms.items()}
        return stms, atos, np.array(eles)

    def eles(self, ctype: str = "mulliken") -> np.ndarray:
        """计算每个原子的电子数

        Args:
            ctype: 电荷类型
        """
        return np.array(self.caler.eles(ctype))

    def eles_obt(self, ctype: str, obt: int) -> np.ndarray:
        """计算指定轨道上的电子数

        Args:
            ctype: 电荷类型
            obt: 分子轨道索引
        """
        return np.array(self.caler.eles_obt(ctype, obt))

    def pi_mocv_v2(
        self, ctype: str = "mulliken"
    ) -> tuple[dict[int, Stm], list[int], np.ndarray]:
        """使用轨道分解方法计算pi电子分布 v2"""
        stms, atos, eles = self.caler.pi_mocv_v2(ctype)
        stms = {atm: Stm(stm) for atm, stm in stms.items()}
        return stms, atos, np.array(eles)


def mulliken(dmat: np.ndarray, smat: np.ndarray, uplo: list[list[int]]) -> np.ndarray:
    """计算Mulliken电荷

    Args:
        dmat: 密度矩阵
        smat: 重叠矩阵
        uplo: 原子轨道索引

    Returns:
        np.ndarray: Mulliken电荷
    """
    return _core.props.atomprop.charge.mulliken(dmat, smat, uplo)  # type: ignore


def lowdin(dmat: np.ndarray, smat: np.ndarray, uplo: list[list[int]]) -> np.ndarray:
    """计算Lowdin电荷

    Args:
        dmat: 密度矩阵
        smat: 重叠矩阵
        uplo: 原子轨道索引

    Returns:
        np.ndarray: Lowdin电荷
    """
    return _core.props.atomprop.charge.lowdin(dmat, smat, uplo)  # type: ignore

def hirshfeld(
    atom_nucs: list[int],
    atom_xyzs: list[list[float]],
    mole_grids: list[list[float]],
    mole_weits: list[float],
    mole_dens: list[float],
    pro_mole_dens: list[float],
) -> list[float]:
    """计算Hirshfeld电荷

    Args:
        atom_nucs: 每个原子的核电荷数
        atom_xyzs: 每个原子的坐标
        mole_grids: 分子格点坐标
        mole_weits: 分子格点权重
        mole_dens: 分子在格点上的密度
        pro_mole_dens: 参考原子密度在格点上的加和

    Returns:
        list[float]: 每个原子的 Hirshfeld 电荷
    """
    return _core.props.atomprop.charge.hirshfeld(  # type: ignore
        atom_nucs, atom_xyzs, mole_grids, mole_weits, mole_dens, pro_mole_dens
    )


def qmat2eles(qmat: np.ndarray, uplo: list[list[int]]) -> list[float]:
    """由原子轨道布居矩阵汇总每个原子的电子数

    Args:
        qmat: 布居矩阵
        uplo: 原子轨道索引

    Returns:
        list[float]: 每个原子的电子数
    """
    return _core.props.atomprop.charge.qmat2eles(qmat, uplo)  # type: ignore