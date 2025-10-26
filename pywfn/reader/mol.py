"""
mol文件读取器
"""

from pywfn import reader
from pywfn.base.geome import Geome
from pywfn.base.basis import Basis
from pywfn.base.coefs import Coefs
from pywfn import core


class MolReader(reader.Reader):
    def __init__(self, path: str) -> None:
        super().__init__(path)
        self.type='mol'
        line3=self.getline(3)
        self.natm=int(line3[:3])
        print(self.natm)
        self.reader=core.reader.MolReader(path) # type: ignore

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
