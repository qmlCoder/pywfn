from pywfn.reader import Reader
from pywfn.base.geome import Geome
from pywfn.base.basis import Basis
from pywfn.base.coefs import Coefs
from pywfn import core

class SdfReader(Reader):
    def __init__(self, path: str, cache: bool = False) -> None:
        super().__init__(path, cache)
        self.type='sdf'
        self.reader=core.reader.SdfReader(path) # type: ignore

    def get_geome(self) -> "Geome":
        geome_core=self.reader.get_geome()
        geome=Geome()
        geome.core=geome_core
        return geome
    
    def get_basis(self)->"Basis":
        basis_core=self.reader.get_basis()
        basis=Basis()
        basis.core=basis_core
        return basis
    
    
    def get_coefs(self)->"Coefs":
        coefs_core=self.reader.get_coefs()
        coefs=Coefs()
        coefs.core=coefs_core
        return coefs

