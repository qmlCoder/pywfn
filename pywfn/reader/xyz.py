"""
本脚本用来读取.xyz文件
xyz文件只包含原子坐标,所以很好读
一个xyz文件可以包含多个分子的坐标，但只生成最后一个坐标的分子对象
"""

from pywfn.base.geome import Geome
from pywfn.base.basis import Basis
from pywfn.base.coefs import Coefs
from pywfn import reader
from pywfn import core

class XyzReader(reader.Reader):
    def __init__(self,path) -> None:
        super().__init__(path)
        self.type='xyz'
        self.path=path
        self.reader=core.reader.XyzReader(path) # type: ignore
    
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
