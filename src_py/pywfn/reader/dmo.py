"""
dmo文件读取器

对应 DMol3 的 .dmo/.dmi 输出文件，除了几何结构外还能读到
原子轨道、键连矩阵与轨道系数。
"""

import numpy as np

from pywfn import _core
from pywfn.base.atoms import Atoms
from pywfn.base.basis import Basis
from pywfn.base.coefs import Coefs


class DmoReader:  # type: ignore
    """DMol3 (.dmo/.dmi) 文件读取器"""

    def __init__(self, _core: _core.reader.DmoReader):  # type: ignore
        """包装 `_core` 的读取器

        Args:
            _core: `_core.reader.DmoReader` 对象
        """
        self.core = _core

    @staticmethod
    def from_path(path: str) -> "DmoReader":
        """由文件路径创建读取器

        Args:
            path: .dmo 文件的路径

        Returns:
            DmoReader: 读取器实例
        """
        return DmoReader(_core.reader.DmoReader(path))  # type: ignore

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

    def read_atms(self) -> tuple[list[int], list[list[float]]]:
        """读取原子序数与坐标

        Returns:
            tuple[list[int], list[list[float]]]: (原子序数列表, 原子坐标列表)
        """
        atms, xyzs = self.core.read_atms()
        return atms, xyzs

    def read_atos(self) -> list[tuple[int, int, str, list[float]]]:
        """读取原子轨道信息

        Returns:
            list[tuple[int, int, str, list[float]]]:
            每个原子轨道的 (原子索引, 壳层索引, 符号, 坐标)
        """
        return self.core.read_atos()

    def read_BM(self) -> np.ndarray:
        """读取键连矩阵

        Returns:
            np.ndarray: 键连矩阵（二维数组）
        """
        return self.core.read_BM()

    def read_coefs(self) -> tuple[list[float], np.ndarray]:
        """读取轨道能量与轨道系数

        Returns:
            tuple[list[float], np.ndarray]: (轨道能量列表, 轨道系数矩阵)
        """
        obt_engs, cmat = self.core.read_coefs()
        return obt_engs, cmat
