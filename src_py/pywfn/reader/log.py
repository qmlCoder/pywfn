"""
此脚本用来提取高斯输出文件的信息
高斯的输出文件包含迭代信息(结构优化、扫描等)
但是封装成分子对象之后就只有一个信息了
所以迭代信息只能在reader对象中存在,且不是默认属性
输出文件中存储的系数可能是球谐的，但是要转为笛卡尔的才方便使用
"""

import numpy as np

from pywfn import _core
from pywfn.base.atoms import Atoms
from pywfn.base.basis import Basis
from pywfn.base.coefs import Coefs


class LogReader:  # type: ignore
    def __init__(self, _core: _core.reader.FchReader):  # type: ignore
        self.core = _core

    @staticmethod
    def from_path(path: str):
        return LogReader(_core.reader.LogReader(path))  # type: ignore

    def get_atoms(self) -> "Atoms":
        return self.core.get_atoms()

    def get_basis(self) -> "Basis":
        return self.core.get_basis()

    def get_coefs(self) -> "Coefs":
        return self.core.get_coefs()

    def read_smat(self) -> "np.ndarray":
        return self.core.read_smat()

    def read_tmat(self) -> "np.ndarray":
        return self.core.read_tmat()

    def read_vmat(self) -> "np.ndarray":
        return self.core.read_vmat()
