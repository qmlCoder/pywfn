"""
本脚本用来读取.xyz文件
xyz文件只包含原子坐标,所以很好读
一个xyz文件可以包含多个分子的坐标，但只生成最后一个坐标的分子对象
"""

from pywfn.base.atoms import Atoms
from pywfn.base.basis import Basis
from pywfn.base.coefs import Coefs
from pywfn import _core


class XyzReader:  # type: ignore
    def __init__(self, _core: _core.reader.XyzReader): # type: ignore
        self.core = _core
    
    @staticmethod
    def from_path(path:str):
        return XyzReader(_core.reader.XyzReader(path)) # type: ignore

    def get_atoms(self) -> "Atoms":
        return self.core.get_atoms()

    def get_basis(self) -> "Basis":
        return self.core.get_basis()

    def get_coefs(self) -> "Coefs":
        return self.core.get_coefs()
