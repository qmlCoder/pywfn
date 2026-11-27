"""
mol文件读取器
"""

from pywfn.base.atoms import Atoms
from pywfn.base.basis import Basis
from pywfn.base.coefs import Coefs
from pywfn import _core


class MolReader():  # type: ignore
    """MDL 分子文件 (.mol) 读取器"""

    def __init__(self, _core: _core.reader.MolReader): # type: ignore
        """包装 `_core` 的读取器

        Args:
            _core: `_core.reader.MolReader` 对象
        """
        self.core = _core
    
    @staticmethod
    def from_path(path: str) -> "MolReader":
        """由文件路径创建读取器

        Args:
            path: .mol 文件的路径

        Returns:
            MolReader: 读取器实例
        """
        return MolReader(_core.reader.MolReader(path)) # type: ignore

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
