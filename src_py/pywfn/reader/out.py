"""
out文件读取器
"""

import numpy as np

from pywfn import _core
from pywfn.base.atoms import Atoms
from pywfn.base.basis import Basis
from pywfn.base.coefs import Coefs


class OutReader:  # type: ignore
    """通用输出文件 (.out) 读取器"""

    def __init__(self, _core: _core.reader.OutReader):  # type: ignore
        """包装 `_core` 的读取器

        Args:
            _core: `_core.reader.OutReader` 对象
        """
        self.core = _core

    @staticmethod
    def from_path(path: str) -> "OutReader":
        """由文件路径创建读取器

        Args:
            path: .out 文件的路径

        Returns:
            OutReader: 读取器实例
        """
        return OutReader(_core.reader.OutReader(path))  # type: ignore

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

    def read_geome(self) -> tuple[list[str], list[list[float]]]:
        """读取几何结构

        Returns:
            tuple[list[str], list[list[float]]]: (原子符号列表, 原子坐标列表)
        """
        syms, xyzs = self.core.read_geome()
        return syms, xyzs

    def read_basis(self) -> dict[int, list[tuple[int, int, float, float]]]:
        """读取基组信息

        Returns:
            dict[int, list[tuple[int, int, float, float]]]:
            {原子索引: [(角动量, 壳层索引, 指数, 系数), ...]}
        """
        return self.core.read_basis()

    def read_coefs(
        self,
    ) -> tuple[list[int], list[int], list[str], list[float], list[float], np.ndarray]:
        """读取分子轨道系数

        Returns:
            tuple[list[int], list[int], list[str], list[float], list[float], np.ndarray]:
            (基函数所属原子索引, 基函数所属壳层索引, 基函数符号,
            轨道能量列表, 轨道占据数列表, 轨道系数矩阵)
        """
        return self.core.read_coefs()

    def read_obtEngs(self) -> tuple[list[float], list[float]]:
        """读取轨道能量和占据数

        Returns:
            tuple[list[float], list[float]]: (轨道能量列表, 轨道占据数列表)
        """
        obt_engs, obt_occs = self.core.read_obtEngs()
        return obt_engs, obt_occs
