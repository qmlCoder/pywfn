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

    def set_cmat(self, form: str, cmat: np.ndarray) -> None:
        """设置系数矩阵

        Args:
            form: 基函数形式, 如 "raw"、"sph"、"car"
            cmat: 轨道系数矩阵
        """
        self.mole.set_cmat(form, cmat)

    def get_cmat(self, form: str = "raw") -> np.ndarray:
        """获取系数矩阵

        统一从 `mole.core` 读取，避免持有克隆的 `core` 时拿到 `set_cmat` 之前的旧值。
        """
        return self.mole.get_cmat(form)

    def get_dmat(self, form: str = "raw") -> np.ndarray:
        """获取密度矩阵（同样从 `mole.core` 读取）"""
        return self.mole.get_dmat(form)

    @property
    def obt_occs(self) -> np.ndarray:
        return self.core.obt_occs

    @property
    def obt_engs(self) -> np.ndarray:
        return self.core.obt_engs

    def obt_idxs(self, if_occ: bool) -> list[int]:
        """获取占据或非占据轨道的索引

        Args:
            if_occ: True 返回占据轨道, False 返回非占据轨道
        """
        return self.core.obt_idxs(if_occ)

    def atm_uplo(self, form: str = "raw") -> np.ndarray:
        return self.mole.basis.atm_uplo(form)

    def is_none(self) -> bool:
        """轨道系数是否为空"""
        return self.core.is_none()

    def obt_eles(self) -> list[int]:
        """占据或虚轨道对应的电子数

        Returns:
            list[int]: [alpha 电子数, beta 电子数]
        """
        return self.core.obt_eles(self.mole.core)

    def equal(self, other: "Coefs") -> bool:
        """判断两个轨道系数是否相同

        Args:
            other: 另一个轨道系数对象
        """
        return self.core.equal(other.core)

    def __repr__(self) -> str:
        return self.core.__repr__()
