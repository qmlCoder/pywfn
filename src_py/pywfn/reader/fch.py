"""
此脚本用来读取fchk文件
fchk文件中有哪些属性是可以用到的？
"""

from pywfn.base.basis import Basis
from pywfn.base.coefs import Coefs
from pywfn.base.atoms import Atoms
from pywfn import _core


class FchReader:  # type: ignore
    def __init__(self, _core: _core.reader.FchReader): # type: ignore
        self.core = _core
    
    @staticmethod
    def from_path(path:str):
        return FchReader(_core.reader.FchReader(path)) # type: ignore

    def get_atoms(self) -> "Atoms":
        return self.core.get_atoms()

    def get_basis(self) -> "Basis":
        return self.core.get_basis()

    def get_coefs(self) -> "Coefs":
        return self.core.get_coefs()
