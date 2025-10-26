"""
定义gjf的读取器，感觉可能也不会用得到
"""
from numpy import ndarray
from pywfn import reader

from pywfn.base.geome import Geome
from pywfn.base.basis import Basis
from pywfn.base.coefs import Coefs

from pywfn.data.elements import elements
from pywfn.data import consts
import re
import numpy as np
from functools import lru_cache
from pywfn import core


class GjfReader(reader.Reader):
    def __init__(self,path) -> None:
        super().__init__(path)
        self.type='gjf'
        self.reader=core.reader.GjfReader(path) # type: ignore

    def get_geome(self) -> Geome:
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
    
    def get_nele(self) -> tuple[int, int]:
        charge,multi=self.read_multi()
        return charge,multi
    
    def get_charge(self) -> int:
        charge,multi=self.read_multi()
        return charge
    
    def get_spin(self) -> int:
        charge,multi=self.read_multi()
        return multi

    @lru_cache
    def read_coord(self):
        finds:list[str]=re.findall(r' ([A-Za-z\d]+) +(-?\d+.\d+) +(-?\d+.\d+) +(-?\d+.\d+)',self.text)
        symbols=[]
        coords=[]
        for s,x,y,z in finds:
            if s.isdigit():
                s=elements[int(s)].symbol
            symbols.append(s)
            coords.append([x,y,z])
        coords=np.array(coords,dtype=np.float32)/consts.Bohr
        return symbols,coords
    
    @lru_cache
    def read_multi(self)->tuple[int,int]:
        """读取电荷与自旋多重度"""
        finds=re.search(r'(-?\d) (\d)',self.text)
        assert finds is not None,'未正确读取电荷与自旋多重度!'
        charge,multi=finds.groups()
        return int(charge),int(multi)