"""
本脚本用来读取.xyz文件
xyz文件只包含原子坐标,所以很好读
一个xyz文件可以包含多个分子的坐标，但只生成最后一个坐标的分子对象
"""

from pywfn.base.atoms import Atoms
from pywfn.base.basis import Basis
from pywfn.base.coefs import Coefs
from pywfn import _core
from pywfn.reader import _wrap_atoms


class XyzReader:  # type: ignore
    """标准 XYZ 坐标文件 (.xyz) 读取器"""

    def __init__(self, _core: _core.reader.XyzReader): # type: ignore
        """包装 `_core` 的读取器

        Args:
            _core: `_core.reader.XyzReader` 对象
        """
        self.core = _core
    
    @staticmethod
    def from_path(path: str) -> "XyzReader":
        """由文件路径创建读取器

        Args:
            path: .xyz 文件的路径

        Returns:
            XyzReader: 读取器实例
        """
        return XyzReader(_core.reader.XyzReader(path)) # type: ignore

    def get_atoms(self) -> "Atoms":
        """获取原子集合（多构象文件里是最后一个构象）

        Returns:
            Atoms: 原子集合（`_core` 层对象）
        """
        return self.core.get_atoms()

    def get_basis(self) -> "Basis | None":
        """获取基组（xyz 文件不含基组，返回 None）

        Returns:
            Basis | None: 恒为 None
        """
        return self.core.get_basis()

    def get_coefs(self) -> "Coefs | None":
        """获取轨道系数（xyz 文件不含系数，返回 None）

        Returns:
            Coefs | None: 恒为 None
        """
        return self.core.get_coefs()

    def read_geomes(self) -> list[Atoms]:
        """读取文件中的所有构象

        每个 `Atoms` 都绑定了一个只含原子的 `Mole`，
        因此可以直接使用 `len()`、`xyzs`、`syms` 等接口。

        Returns:
            list[Atoms]: 每个构象对应一个原子集合
        """
        return [_wrap_atoms(atoms) for atoms in self.core.read_geomes()]
