"""
读取.molden文件的脚本
"""
from pywfn.base import Atoms, Basis, Coefs
from pywfn import _core


sym2ang = {'s': 0, 'p': 1, 'd': 2, 'f': 3, 'g': 4, 'h': 5, 'i': 6}


class MdeReader:  # type: ignore
    def __init__(self, _core: _core.reader.MdeReader): # type: ignore
        self.core = _core

    @staticmethod
    def from_path(path:str):
        return MdeReader(_core.reader.MdeReader(path)) # type: ignore

    def get_atoms(self) -> "Atoms":
        return self.core.get_atoms()

    def get_basis(self) -> "Basis":
        return self.core.get_basis()

    def get_coefs(self) -> "Coefs":
        return self.core.get_coefs()
