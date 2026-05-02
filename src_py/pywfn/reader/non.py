"""
不从文件中读取信息的读取器
"""
import numpy as np

from pywfn.base.atoms import Atoms
from pywfn.base.basis import Basis
from pywfn.base.coefs import Coefs

from pywfn import _core


class NonReader:  # type: ignore
    def __init__(self, _core: _core.reader.NonReader): # type: ignore
        self.core = _core
    
    @staticmethod
    def from_path(path:str):
        return NonReader(_core.reader.NonReader(path)) # type: ignore

    def get_atoms(self) -> "Atoms":
        return self.core.get_atoms()

    def get_basis(self) -> "Basis":
        return self.core.get_basis()

    def get_coefs(self) -> "Coefs":
        return self.core.get_coefs()
