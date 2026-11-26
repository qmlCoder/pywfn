"""
out文件读取器
"""

from pywfn import _core
from pywfn.base.atoms import Atoms
from pywfn.base.basis import Basis
from pywfn.base.coefs import Coefs


class OutReader:  # type: ignore
    def __init__(self, _core: _core.reader.OutReader):  # type: ignore
        self.core = _core

    @staticmethod
    def from_path(path: str):
        return OutReader(_core.reader.OutReader(path))  # type: ignore

    def get_atoms(self) -> "Atoms":
        return self.core.get_atoms()

    def get_basis(self) -> "Basis | None":
        return self.core.get_basis()

    def get_coefs(self) -> "Coefs | None":
        return self.core.get_coefs()
