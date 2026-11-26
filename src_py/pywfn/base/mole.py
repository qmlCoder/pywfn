"""
基础的分子对象
"""

import numpy as np

from pywfn import _core
from pywfn.base.atoms import Atoms
from pywfn.base.basis import Basis
from pywfn.base.coefs import Coefs


class Mole:  # type: ignore
    """基础的分子对象"""

    def __init__(self, _core: _core.base.Mole) -> None:  # type: ignore
        self.core = _core

    @staticmethod
    def from_file(path: str) -> "Mole":
        """从文件读取分子结构"""
        return Mole(_core.base.Mole(path))  # type: ignore

    @property
    def atoms(self) -> Atoms:
        """原子集合（每次访问都从 core 重新获取，保证与 `mole.core` 一致）"""
        return Atoms(self.core.atoms(), self)

    @property
    def basis(self) -> Basis:
        """基组信息（每次访问都从 core 重新获取）"""
        return Basis(self.core.basis(), self)

    @property
    def coefs(self) -> Coefs:
        """轨道系数（每次访问都从 core 重新获取，`set_cmat` 后能读到新值）"""
        return Coefs(self.core.coefs(), self)

    @property
    def charge(self) -> int:
        return self.core.charge()

    @property
    def multip(self) -> int:
        return self.core.multip()

    def syms(self) -> list[str]:
        """所有原子的元素符号列表"""
        return self.core.syms()

    def xyzs(self) -> np.ndarray:
        """所有原子的笛卡尔坐标列表"""
        return np.array(self.core.xyzs())

    def __repr__(self) -> str:
        return self.core.__repr__()

    def get_cmat(self, form: str) -> np.ndarray:
        """获取系数矩阵"""
        return self.core.get_cmat(form)

    def get_smat(self, form: str) -> np.ndarray:
        """获取重叠矩阵"""
        return self.core.get_smat(form)

    def get_dmat(self, form: str) -> np.ndarray:
        """获取密度矩阵"""
        return self.core.get_dmat(form)

    def set_cmat(self, form: str, cmat: np.ndarray) -> None:
        """设置系数矩阵

        Args:
            form: 基函数形式, 如 "raw"、"sph"、"car"
            cmat: 轨道系数矩阵
        """
        self.core.set_cmat(cmat, form)

    def border(self, bord: float) -> np.ndarray:
        """获取分子的边界框"""
        return self.core.border(bord)

    def is_open(self) -> bool:
        """判断分子是否为开壳层"""
        return self.core.is_open()

    def clone(self) -> 'Mole':
        """克隆分子对象"""
        return Mole(self.core.clone())
