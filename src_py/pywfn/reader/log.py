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
from pywfn.reader import GauJob, Nics, Vibe
from pywfn.reader import _wrap_atoms, _wrap_gau_job, _wrap_vibe


class LogReader:  # type: ignore
    """Gaussian 输出文件 (.log/.out) 读取器"""

    def __init__(self, _core: _core.reader.LogReader):  # type: ignore
        """包装 `_core` 的读取器

        Args:
            _core: `_core.reader.LogReader` 对象
        """
        self.core = _core

    @staticmethod
    def from_path(path: str) -> "LogReader":
        """由文件路径创建读取器

        Args:
            path: .log/.out 文件的路径

        Returns:
            LogReader: 读取器实例
        """
        return LogReader(_core.reader.LogReader(path))  # type: ignore

    def get_atoms(self) -> "Atoms":
        """获取原子集合

        Returns:
            Atoms: 原子集合（`_core` 层对象）
        """
        return self.core.get_atoms()

    def get_basis(self) -> "Basis | None":
        """获取基组

        Returns:
            Basis | None: 基组（`_core` 层对象），文件中没有时为 None
        """
        return self.core.get_basis()

    def get_coefs(self) -> "Coefs | None":
        """获取轨道系数

        Returns:
            Coefs | None: 轨道系数（`_core` 层对象），文件中没有时为 None
        """
        return self.core.get_coefs()

    def read_job(self) -> str:
        """从输出文件中读取任务信息"""
        return self.core.read_job()

    def read_gau_job(self) -> GauJob | None:
        """解析作业关键词，得到任务类型、泛函与基组

        Returns:
            GauJob | None: 解析成功时返回任务信息，否则返回 None
        """
        core_job = self.core.read_gau_job()
        return None if core_job is None else _wrap_gau_job(core_job)

    def _read_coefs(
        self, spin: str
    ) -> tuple[
        list[int], list[int], list[str], list[float], list[str], np.ndarray
    ]:
        """读取指定自旋的系数矩阵

        Args:
            spin: 自旋类型，如 "Alpha" 或 "Beta"

        Returns:
            tuple[list[int], list[int], list[str], list[float], list[str], np.ndarray]:
            (基函数所属原子索引, 基函数所属壳层索引, 基函数符号,
            轨道能量列表, 轨道相关的字符串列表, 轨道系数矩阵)
        """
        return self.core._read_coefs(spin)

    def read_coefs(
        self,
    ) -> tuple[
        list[int], list[int], list[str], list[float], list[str], np.ndarray
    ] | None:
        """读取轨道系数（结果会被缓存）

        Returns:
            tuple | None: (基函数所属原子索引, 基函数所属壳层索引, 基函数符号,
            轨道能量列表, 轨道相关的字符串列表, 轨道系数矩阵)，文件中没有时返回 None
        """
        return self.core.read_coefs()

    def read_geome(
        self, irow: int | None = None
    ) -> tuple[list[int], list[list[float]]]:
        """读取分子几何结构

        Args:
            irow: 起始行号，为 None 时自动查找

        Returns:
            tuple[list[int], list[list[float]]]: (原子序数列表, 原子坐标列表)
        """
        atms, xyzs = self.core.read_geome(irow)
        return atms, xyzs

    def read_geomes(self) -> list[Atoms]:
        """读取所有优化步的几何结构

        Returns:
            list[Atoms]: 每一步结构对应一个原子集合
        """
        return [_wrap_atoms(atoms) for atoms in self.core.read_geomes()]

    def read_energys(self) -> list[float] | None:
        """读取热力学能量数据

        Returns:
            list[float] | None: 能量列表，没有能量信息时返回 None
        """
        return self.core.read_energys()

    def read_freqs(self) -> list[float] | None:
        """读取振动频率

        Returns:
            list[float] | None: 频率列表，没有频率信息时返回 None
        """
        return self.core.read_freqs()

    def read_symMat(self, row: int, loc: int) -> np.ndarray:
        """读取对称矩阵

        高斯输出里的矩阵是下三角存储的，需要给出起始行列才能正确解析。

        Args:
            row: 读取数据的起始行号
            loc: 读取数据的起始列号

        Returns:
            np.ndarray: 二维矩阵
        """
        return self.core.read_symMat(row, loc)

    def read_dmat(self) -> np.ndarray | None:
        """读取密度矩阵

        Returns:
            np.ndarray | None: 密度矩阵，不存在时返回 None
        """
        return self.core.read_dmat()

    def read_smat(self) -> "np.ndarray | None":
        """从输出文件中读取重叠矩阵"""
        return self.core.read_smat()

    def read_tmat(self) -> "np.ndarray | None":
        """从输出文件中读取动能矩阵"""
        return self.core.read_tmat()

    def read_vmat(self) -> "np.ndarray | None":
        """从输出文件中读取势能矩阵"""
        return self.core.read_vmat()

    def read_ecps(self) -> dict[int, int] | None:
        """读取赝势基组的有效核电荷数

        Returns:
            dict[int, int] | None: {原子索引: 有效核电荷数}，不存在时返回 None
        """
        return self.core.read_ecps()

    def read_vibrations(self) -> list[Vibe] | None:
        """读取振动模式信息

        Returns:
            list[Vibe] | None: 振动模式列表，文件中不含振动信息时返回 None
        """
        core_vibes = self.core.read_vibrations()
        return None if core_vibes is None else [_wrap_vibe(v) for v in core_vibes]

    def read_nics(self) -> list["Nics"]|None:
        """从输出文件中读取NICS数据"""
        cores = self.core.read_nics()
        return None if cores is None else [Nics(core) for core in cores]

    def read_normalend(self) -> bool:
        """判断输出文件是否正常结束

        Returns:
            bool: 正常结束时返回 True，否则返回 False
        """
        return self.core.read_normalend()
