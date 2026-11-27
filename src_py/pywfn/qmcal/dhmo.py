import numpy as np

from pywfn import _core
from pywfn.base.mole import Mole


class DmoRes:
    """dHMO 计算结果"""

    def __init__(self, core) -> None:  # type: ignore[no-untyped-def]
        self._core = core

    @property
    def atos(self) -> list[tuple[int, int, str, list[float]]]:
        """原子轨道列表，每个元素为 (原子索引, 壳层索引, 轨道符号, 方向向量 [x,y,z])"""
        return self._core.atos

    @property
    def bm(self) -> np.ndarray:
        """键连矩阵"""
        return self._core.bm()

    @property
    def cm(self) -> np.ndarray:
        """轨道系数矩阵"""
        return self._core.cm()

    @property
    def BM(self) -> np.ndarray:
        """键连矩阵（对应上游字段 `BM`，与 `bm` 内容一致）"""
        return self._core.BM

    @property
    def CM(self) -> np.ndarray:
        """轨道系数矩阵（对应上游字段 `CM`，与 `cm` 内容一致）"""
        return self._core.CM

    @property
    def es(self) -> list[float]:
        """分子轨道能量"""
        return self._core.es


def effective_charge(mole: Mole, atm: int) -> float:
    """计算 sp2 杂化原子的有效核电荷数"""
    return _core.qmcal.dhmo.effective_charge(mole.core, atm)


def overlap_integral(mole: Mole, ai: int, aj: int) -> float:
    """计算两个原子间主量子数最高的 Slater 型 p 轨道的重叠积分"""
    return _core.qmcal.dhmo.overlap_integral(mole.core, ai, aj)


def free_valence(mole: Mole, atm: int, sum_order: float) -> float:
    """计算原子的自由价"""
    return _core.qmcal.dhmo.free_valence(mole.core, atm, sum_order)


def dhmo_calc(path: str) -> tuple[Mole, DmoRes]:
    """执行 dHMO 计算，返回分子结构与结果"""
    mole, res = _core.qmcal.dhmo.dhmo_calc(path)
    return Mole(mole), DmoRes(res)


def output(dmi_path: str, dmo_path: str | None = None) -> None:
    """计算并输出 dHMO 结果到 .dmo 文件"""
    _core.qmcal.dhmo.output(dmi_path, dmo_path)


def HMO_Analysis(
    mole: Mole,
    bm: np.ndarray,
    cm: np.ndarray,
    ato_atms: list[int],
    nele: float,
) -> tuple[np.ndarray, list[float], float]:
    """由 Hückel 哈密顿矩阵求 dHMO 的轨道能量与相关量

    Args:
        mole: 分子对象
        bm: 键连矩阵，形状为 (nato, nato)
        cm: 轨道系数矩阵，形状为 (nato, nato)
        ato_atms: 每个原子轨道所属的原子索引
        nele: 电子数

    Returns:
        tuple[np.ndarray, list[float], float]: (密度矩阵, 轨道能量列表, 电子能量)
    """
    return _core.qmcal.dhmo.HMO_Analysis(mole.core, bm, cm, ato_atms, nele)  # type: ignore
