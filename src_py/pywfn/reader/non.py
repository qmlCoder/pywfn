"""
不从文件中读取信息的读取器
"""

from pywfn import _core
from pywfn.base.atoms import Atoms
from pywfn.base.basis import Basis
from pywfn.base.coefs import Coefs


class NonReader:  # type: ignore
    def __init__(self, _core: _core.reader.NonReader):  # type: ignore
        self.core = _core

    @staticmethod
    def from_path(path: str):
        return NonReader(_core.reader.NonReader(path))  # type: ignore

    def set_atoms(self, atms: list[int], xyzs: list[list[float]]) -> None:
        """设置原子索引和坐标

        Args:
            atms: 原子索引列表
            xyzs: 原子坐标列表，每个元素为 [x, y, z] 三元组
        """
        self.core.set_atoms(atms, xyzs)

    def get_atoms(self) -> "Atoms":
        """获取原子集合（与其它读取器保持一致的接口）"""
        return self.core.get_atoms()

    def get_basis(self) -> "Basis | None":
        """获取基组信息（NonReader 无基组，返回 None）"""
        return self.core.get_basis()

    def get_coefs(self) -> "Coefs | None":
        """获取轨道系数（NonReader 无系数，返回 None）"""
        return self.core.get_coefs()
