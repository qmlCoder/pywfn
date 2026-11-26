from pywfn import _core
from pywfn.base.atoms import Atoms
from pywfn.base.basis import Basis
from pywfn.base.coefs import Coefs


class CubReader:  # type: ignore
    def __init__(self, _core: _core.reader.CubReader):  # type: ignore
        self.core = _core

    @staticmethod
    def from_path(path: str):
        return CubReader(_core.reader.CubReader(path))  # type: ignore

    def get_atoms(self) -> "Atoms":
        """获取分子几何信息"""
        return self.core.get_atoms()

    def get_basis(self) -> "Basis | None":
        """获取基组信息 (cub文件可能不包含基组信息)"""
        return self.core.get_basis()

    def get_coefs(self) -> "Coefs | None":
        """获取轨道系数 (cub文件可能不包含轨道系数)"""
        return self.core.get_coefs()
