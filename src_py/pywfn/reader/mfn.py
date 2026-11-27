"""
mfn/mwfn 文件读取器

mfn（Multiwfn 波函数）文件只保存波函数本身：原子、基函数、轨道系数，
不包含成键信息与格点数据。这里的方法与 `_core.reader.MfnReader` 一一对应。

说明：`get_atoms`/`get_basis`/`get_coefs` 和其它 Reader 保持一致，
直接返回 `_core` 层对象；需要绑定分子的 Python 包装对象时请用 `Mole.from_file`。
"""

import numpy as np

from pywfn import _core
from pywfn.reader import BasData, _wrap_bas_data


class MfnReader:  # type: ignore
    """Multiwfn 波函数文件 (.mwfn) 读取器"""

    def __init__(self, _core: "_core.reader.MfnReader") -> None:  # type: ignore
        """包装 `_core` 的读取器

        Args:
            _core: `_core.reader.MfnReader` 对象
        """
        self.core = _core

    @staticmethod
    def from_path(path: str) -> "MfnReader":
        """由文件路径创建读取器

        Args:
            path: .mwfn 文件的路径

        Returns:
            MfnReader: 读取器实例
        """
        return MfnReader(_core.reader.MfnReader(path))  # type: ignore

    def get_atoms(self) -> "_core.base.Atoms":
        """获取原子集合

        Returns:
            _core.base.Atoms: 原子集合（`_core` 层对象）
        """
        return self.core.get_atoms()

    def get_basis(self) -> "_core.base.Basis | None":
        """获取基组

        Returns:
            _core.base.Basis | None: 基组（`_core` 层对象），文件不含基组时为 None
        """
        return self.core.get_basis()

    def get_coefs(self) -> "_core.base.Coefs | None":
        """获取轨道系数

        Returns:
            _core.base.Coefs | None: 轨道系数（`_core` 层对象），文件不含系数时为 None
        """
        return self.core.get_coefs()

    def read_geome(self) -> tuple[list[int], list[list[float]], list[float]]:
        """读取几何结构

        内部坐标使用 Bohr，有效核电荷数用于标记赝势。

        Returns:
            tuple[list[int], list[list[float]], list[float]]:
            (原子序数列表, 原子坐标列表, 有效核电荷数列表)
        """
        return self.core.read_geome()

    def read_basis(
        self,
    ) -> tuple[list[int], list[int], list[str], list[BasData], list[str]]:
        """读取基函数信息

        Returns:
            tuple[list[int], list[int], list[str], list[BasData], list[str]]:
            (基函数所属原子索引, 基函数所属壳层索引, 基函数符号,
            基函数原始数据列表, 壳层基函数形式列表)
        """
        ato_atms, ato_shls, ato_syms, ato_bass, ato_forms = self.core.read_basis()
        datas = [_wrap_bas_data(bas) for bas in ato_bass]
        return ato_atms, ato_shls, ato_syms, datas, ato_forms

    def read_coefs(
        self,
    ) -> tuple[
        list[int], list[int], list[str], list[float], list[float], np.ndarray
    ]:
        """读取轨道系数、能量和占据数

        Returns:
            tuple[list[int], list[int], list[str], list[float], list[float], np.ndarray]:
            (基函数所属原子索引, 基函数所属壳层索引, 基函数符号,
            轨道能量列表, 轨道占据数列表, 轨道系数矩阵)
        """
        return self.core.read_coefs()

    def read_neles(self) -> tuple[int, int]:
        """读取 α 和 β 电子数

        Returns:
            tuple[int, int]: (α 电子数, β 电子数)
        """
        nela, nelb = self.core.read_neles()
        return nela, nelb
