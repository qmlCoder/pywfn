"""
计算分子轨道空间波函数
基函数的波函数
原子轨道的波函数
原子的波函数
"""
from pywfn.base import Mole
import numpy as np
from pywfn import _core

Array = np.ndarray


class Calculator:
    def __init__(self, mole: Mole) -> None:
        self.mole = mole
        self.caler = _core.gridprop.wfnfunc.Calculator(mole.core)  # type: ignore # 核心计算器

    def obt_wfn(self, obt: int, grids: np.ndarray, level: int):
        val0, val1, val2 = self.caler.obt_wfn(obt, grids, level)
        val0 = np.array(val0)
        val1 = np.array(val1)
        val2 = np.array(val2)
        return val0, val1, val2

    # 所有原子轨道的波函数
    def ato_wfn(self, grids: np.ndarray, level: int, atms: list[int] | None = None):
        """计算所有原子轨道"""
        return self.caler.ato_wfn(grids, level, atms)

    def fragOverlap(self, mols: list[Mole], obts: list[int]):
        """
        计算一个分子两个片段之间的重叠
        将一个分子的两个片段单独导出并进行单点计算
        指定这两个分子的分子轨道索引
        """
        from pywfn.gridprop import dftgrid
        grids, weits = dftgrid.Calculator(self.mole).mol_grids()
        mol0, mol1 = mols
        obt0, obt1 = obts
        wfn0 = Calculator(mol0).obt_wfn(obt0, grids, 0)[0]
        wfn1 = Calculator(mol1).obt_wfn(obt1, grids, 0)[0]
        return np.sum(wfn0*wfn1*weits)


def gtf(xyzs: np.ndarray, lmn: list[int], alp: float, level: int):
    vals = _core.gridprop.wfnfunc.gtf(xyzs, lmn, alp, level)  # type: ignore
    return vals
