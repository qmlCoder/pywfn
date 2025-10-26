"""
不从文件中读取信息的读取器
"""
import numpy as np

from pywfn import reader

from pywfn.base.geome import Geome
from pywfn.base.basis import Basis
from pywfn.base.coefs import Coefs

from pywfn import core

class NonReader(reader.Reader):
    def __init__(self,path:str=""):
        super().__init__(path)
        self.type:str="any"
        self.charge:int=0
        self.spin:int=1
        self.geome=Geome()
        self.basis=Basis()
        self.coefs=Coefs()
        self.reader=core.reader.NonReader(path) # type: ignore

    def get_charge(self) -> int:
        return self.charge
    
    def get_spin(self) -> int:
        return self.spin

    def get_geome(self) -> Geome:
        return self.geome

    def get_basis(self) -> Basis:
        return self.basis

    def get_coefs(self) -> Coefs:
        return self.coefs