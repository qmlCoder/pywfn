"""
基础的分子对象，其属性应该是标准的，已知的
一个分子对象应该有哪些属性？基本属性(必须属性),计算属性(需要计算才能的到的属性)
结构相关
- 原子
- 键
- 法向量
- 重原子

分子轨道相关
轨道数量
轨道类型
是否分为α和β
占据轨道序数
非占据轨道序数
重叠矩阵

读取器
"""

from functools import cached_property
from pywfn.base.coefs import Coefs
from pywfn.base.basis import Basis
from pywfn.base.atoms import Atoms
from pywfn import _core

import numpy as np


class Mole:  # type: ignore
    """基础的分子对象"""

    def __init__(self, _core: _core.base.Mole) -> None:  # type: ignore
        self.core = _core

    @staticmethod
    def from_file(path: str) -> "Mole":
        return Mole(_core.base.Mole(path))  # type: ignore

    @cached_property
    def atoms(self) -> Atoms:
        return Atoms(self.core.atoms())

    @cached_property
    def basis(self) -> Basis:
        return Basis(self.core.basis())

    @cached_property
    def coefs(self) -> Coefs:
        return Coefs(self.core.coefs())

    def __repr__(self) -> str:
        return self.core.__repr__()

    def get_cmat(self, form: str) -> np.ndarray:
        return self.core.get_cmat(form)

    def set_cmat(self, form: str, cmat):
        self.core.set_cmat(form, cmat)

    def border(self):
        return self.core.border()
