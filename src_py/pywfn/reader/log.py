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
from pywfn.reader import Nics


class LogReader:  # type: ignore

    def __init__(self, _core: _core.reader.LogReader):  # type: ignore
        self.core = _core

    @staticmethod
    def from_path(path: str):
        return LogReader(_core.reader.LogReader(path))  # type: ignore

    def get_atoms(self) -> "Atoms":
        return self.core.get_atoms()

    def get_basis(self) -> "Basis | None":
        return self.core.get_basis()

    def get_coefs(self) -> "Coefs | None":
        return self.core.get_coefs()

    def read_job(self) -> str:
        """从输出文件中读取任务信息"""
        return self.core.read_job()

    def read_smat(self) -> "np.ndarray | None":
        """从输出文件中读取重叠矩阵"""
        return self.core.read_smat()

    def read_tmat(self) -> "np.ndarray | None":
        """从输出文件中读取动能矩阵"""
        return self.core.read_tmat()

    def read_vmat(self) -> "np.ndarray | None":
        """从输出文件中读取势能矩阵"""
        return self.core.read_vmat()

    def read_nics(self) -> list["Nics"]|None:
        """从输出文件中读取NICS数据"""
        cores = self.core.read_nics()
        return None if cores is None else [Nics(core) for core in cores]
