from pywfn.base.atoms import Atoms
from pywfn.base.basis import Basis
from pywfn.base.coefs import Coefs
from pywfn import _core


class SdfReader:  # type: ignore
    def __init__(self, _core: _core.reader.SdfReader): # type: ignore
        self.core = _core
    
    @staticmethod
    def from_path(path:str):
        return SdfReader(_core.reader.SdfReader(path)) # type: ignore

    def get_atoms(self) -> "Atoms":
        return self.core.get_atoms()

    def get_basis(self) -> "Basis":
        return self.core.get_basis()

    def get_coefs(self) -> "Coefs":
        return self.core.get_coefs()
