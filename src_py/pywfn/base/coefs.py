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
from pywfn import _core


class Coefs:

    def __init__(self, core: _core.base.Coefs) -> None:  # type: ignore
        self.core = core

    @staticmethod
    def new(ato_atms, ato_shls, ato_syms, obt_engs, obt_occs, cmat) -> "Coefs":
        core = _core.base.Coefs(ato_atms, ato_shls, ato_syms, obt_engs, obt_occs, cmat)  # type: ignore
        return Coefs(core)

    def get_cmat(self, form: str) -> np.ndarray:
        return np.array(self.core.get_cmat(form))

    def __repr__(self) -> str:
        return self.core.__repr__()
