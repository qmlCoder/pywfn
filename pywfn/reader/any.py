"""
不从文件中读取信息的读取器
"""
import numpy as np

from pywfn import reader
from pywfn.base.basis import Basis
from pywfn import base
class AnyReader(reader.Reader):
    def __init__(self,path:str=''):
        super().__init__(path)
        self.type:str='any'
        self.charge:int=0
        self.spin:int=1
        self.geome:"base.Geome|None"=None
        self.basis:"base.Basis|None"=None
        self.coefs:"base.Coefs|None"=None
        
    
    def get_charge(self) -> int:
        return self.charge
    
    def get_spin(self) -> int:
        return self.spin

    def get_geome(self) -> "base.Geome":
        assert self.geome is not None,"未设置结构信息"
        return self.geome
    
    def get_basis(self) -> "base.Basis":
        assert self.basis is not None,"未设置机组信息"
        return self.basis
    
    def get_coefs(self) -> "base.Coefs":
        assert self.coefs is not None,"未设置轨道信息"
        return self.coefs