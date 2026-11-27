"""
inp文件读取器
"""

from pywfn import _core
from pywfn.base.atoms import Atoms
from pywfn.base.basis import Basis
from pywfn.base.coefs import Coefs


class InpReader:  # type: ignore
    """通用输入文件 (.inp) 读取器"""

    def __init__(self, _core: _core.reader.InpReader):  # type: ignore
        """包装 `_core` 的读取器

        Args:
            _core: `_core.reader.InpReader` 对象
        """
        self.core = _core

    @staticmethod
    def from_path(path: str) -> "InpReader":
        """由文件路径创建读取器

        Args:
            path: .inp 文件的路径

        Returns:
            InpReader: 读取器实例
        """
        return InpReader(_core.reader.InpReader(path))  # type: ignore

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

        与其它读取器不同，这里的结构是从输入文件里直接抄的，
        所以拿到的是元素符号而不是原子序数。

        Returns:
            tuple[list[str], list[list[float]]]: (原子符号列表, 原子坐标列表)
        """
        syms, xyzs = self.core.read_geome()
        return syms, xyzs
