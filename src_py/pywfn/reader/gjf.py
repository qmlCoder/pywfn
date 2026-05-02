"""
定义gjf的读取器，感觉可能也不会用得到
"""
from numpy import ndarray

from pywfn.base.atoms import Atoms
from pywfn.base.basis import Basis
from pywfn.base.coefs import Coefs
from pywfn import _core


class GjfReader:  # type: ignore
    def __init__(self, _core: _core.reader.GjfReader): # type: ignore
        self.core = _core
    
    @staticmethod
    def from_path(path:str):
        return GjfReader(_core.reader.GjfReader(path)) # type: ignore

    def get_atoms(self) -> "Atoms":
        return self.core.get_atoms()

    def get_basis(self) -> "Basis":
        return self.core.get_basis()

    def get_coefs(self) -> "Coefs":
        return self.core.get_coefs()
