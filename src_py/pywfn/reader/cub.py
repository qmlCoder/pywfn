from pywfn.base.basis import Basis
from pywfn.base.coefs import Coefs
from pywfn.base.atoms import Atoms
from pywfn.base import Mole
from pywfn import _core


class CubReader:  # type: ignore
    def __init__(self, _core: _core.reader.CubReader): # type: ignore
        self.core = _core
    
    @staticmethod
    def from_path(path:str):
        return CubReader(_core.reader.CubReader(path)) # type: ignore

    def get_atoms(self) -> "Atoms":
        """获取分子几何信息"""
        print(self.core)
        return Atoms(self.core.get_atoms())

    def get_basis(self) -> None:
        return None

    def get_coefs(self) -> None:
        return None
