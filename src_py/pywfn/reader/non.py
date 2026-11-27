"""
不从文件中读取信息的读取器
"""

from pywfn import _core
from pywfn.base.atoms import Atoms
from pywfn.base.basis import Basis
from pywfn.base.coefs import Coefs


class NonReader:  # type: ignore
    """空读取器：用于手动构建原子信息，不从文件中读取"""

    def __init__(self, _core: _core.reader.NonReader):  # type: ignore
        """包装 `_core` 的读取器

        Args:
            _core: `_core.reader.NonReader` 对象
        """
        self.core = _core
        # `_core` 的 NonReader 只暴露 set_atoms，没有取回原子的接口，
        # 因此在 Python 侧同时记一份，供 get_atoms 重建原子集合使用
        self._atms: list[int] = []
        self._xyzs: list[list[float]] = []

    @staticmethod
    def from_path(path: str) -> "NonReader":
        """创建一个空的读取器

        Args:
            path: 路径（仅用于占位，不会读取文件内容）

        Returns:
            NonReader: 读取器实例
        """
        return NonReader(_core.reader.NonReader(path))  # type: ignore

    def set_atoms(self, atms: list[int], xyzs: list[list[float]]) -> None:
        """设置原子索引和坐标

        Args:
            atms: 原子索引列表（同时作为核电荷数使用）
            xyzs: 原子坐标列表，每个元素为 [x, y, z] 三元组
        """
        self.core.set_atoms(atms, xyzs)
        self._atms = list(atms)
        self._xyzs = [list(xyz) for xyz in xyzs]

    def get_atoms(self) -> "_core.base.Atoms":
        """获取原子集合（与其它读取器保持一致的接口）

        `_core` 的 NonReader 没有 get_atoms，这里用 `set_atoms` 记录的数据重建，
        重建规则与 Rust 侧一致：核电荷数取原子序数。

        Returns:
            _core.base.Atoms: 原子集合（`_core` 层对象）
        """
        crgs = [float(atm) for atm in self._atms]
        return _core.base.Atoms(self._atms, self._xyzs, crgs)  # type: ignore

    def get_basis(self) -> "Basis | None":
        """获取基组信息（NonReader 无基组，返回 None）

        Returns:
            Basis | None: 恒为 None
        """
        return None

    def get_coefs(self) -> "Coefs | None":
        """获取轨道系数（NonReader 无系数，返回 None）

        Returns:
            Coefs | None: 恒为 None
        """
        return None
