"""
存储分子轨道系数矩阵
存储直接读到的系数矩阵及相关数据
- CM 系数矩阵
- atoAtms 每一行对应的原子
- atoShls 每一个原子轨道的角动量
- atoSyms 每个原子轨道的符号
- obtEngs 分子轨道能量
- obtOccs 分子轨道是否占据

轨道的类型要能够在混合、全笛卡尔和全球谐之间转换
"""

import numpy as np

import pywfn
from pywfn import _core


class Coefs:
    def __init__(self, core: _core.base.Coefs, mole: "pywfn.base.Mole") -> None:  # type: ignore
        self.core = core
        self.mole = mole

    def set_cmat(self, cmat: np.ndarray, form: str = "raw") -> None:
        basis = self.mole.basis
        self.core.set_cmat(form, cmat, basis.core)
        self.mole.set_cmat(cmat, form)

    def get_cmat(self, form: str = "raw") -> np.ndarray:
        basis = self.mole.basis
        return self.core.get_cmat(form, basis.core)

    def get_dmat(self, form: str = "raw") -> np.ndarray:
        basis = self.mole.basis
        return self.core.get_dmat(form, basis.core)

    @property
    def obt_occs(self) -> np.ndarray:
        return self.core.obt_occs

    @property
    def obt_engs(self) -> np.ndarray:
        return self.core.obt_engs

    def obt_idxs(self, if_occ: bool) -> list[int]:
        return self.core.obt_idxs(if_occ)

    def atm_uplo(self, form: str = "raw") -> np.ndarray:
        return self.core.atm_uplo(form)

    def __repr__(self) -> str:
        return self.core.__repr__()
