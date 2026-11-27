"""
读取.molden文件的脚本
"""
from pywfn.base import Atoms, Basis, Coefs
from pywfn import _core
from pywfn.reader import BasData, _wrap_bas_data


sym2ang = {'s': 0, 'p': 1, 'd': 2, 'f': 3, 'g': 4, 'h': 5, 'i': 6}


class MdeReader:  # type: ignore
    """Molden 波函数文件 (.molden) 读取器"""

    def __init__(self, _core: _core.reader.MdeReader): # type: ignore
        """包装 `_core` 的读取器

        Args:
            _core: `_core.reader.MdeReader` 对象
        """
        self.core = _core

    @staticmethod
    def from_path(path: str) -> "MdeReader":
        """由文件路径创建读取器

        Args:
            path: .molden 文件的路径

        Returns:
            MdeReader: 读取器实例
        """
        return MdeReader(_core.reader.MdeReader(path)) # type: ignore

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

    def read_geome(self) -> tuple[list[int], list[list[float]]]:
        """读取几何结构

        Returns:
            tuple[list[int], list[list[float]]]: (原子序数列表, 原子坐标列表)
        """
        atms, xyzs = self.core.read_geome()
        return atms, xyzs

    def read_basis(self) -> tuple[list[int], list[int], list[str], list[BasData]]:
        """读取基函数信息

        Returns:
            tuple[list[int], list[int], list[str], list[BasData]]:
            (基函数所属原子索引, 基函数所属壳层索引, 基函数符号,
            基函数原始数据列表)
        """
        ato_atms, ato_shls, ato_syms, ato_bass = self.core.read_basis()
        datas = [_wrap_bas_data(bas) for bas in ato_bass]
        return ato_atms, ato_shls, ato_syms, datas

    def read_coefs(self) -> tuple[list[float], list[float], list[list[float]]]:
        """读取轨道系数

        Returns:
            tuple[list[float], list[float], list[list[float]]]:
            (轨道能量列表, 轨道占据数列表, 轨道系数矩阵)
        """
        obt_engs, obt_occs, cmat = self.core.read_coefs()
        return obt_engs, obt_occs, cmat

    def read_neles(self) -> tuple[int, int]:
        """读取电子数

        Returns:
            tuple[int, int]: (α 电子数, β 电子数)
        """
        nela, nelb = self.core.read_neles()
        return nela, nelb
