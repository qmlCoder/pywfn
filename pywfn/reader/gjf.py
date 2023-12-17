"""
定义gjf的读取器，感觉可能也不会用得到
"""
from numpy import ndarray
from pywfn.reader import lutils
import re
import numpy as np
from functools import lru_cache

class GjfReader(lutils.Reader):
    def __init__(self,path) -> None:
        super().__init__(path)
    
    def get_symbols(self) -> list[str]:
        symbols,coords=self.read_coord()
        return symbols
    
    def get_coords(self) -> ndarray:
        symbols,coords=self.read_coord()
        return coords
    
    def get_charge(self) -> int:
        charge,multi=self.read_multi()
        return charge
    
    def get_spin(self) -> int:
        charge,multi=self.read_multi()
        return multi

    @lru_cache
    def read_coord(self):
        finds=re.findall(' ([A-Za-z]+) +(-?\d+.\d+) +(-?\d+.\d+) +(-?\d+.\d+)',self.text)
        symbols=[]
        coords=[]
        for s,x,y,z in finds:
            symbols.append(s)
            coords.append([x,y,z])
        return symbols,np.array(coords,dtype=float)
    
    @lru_cache
    def read_multi(self)->tuple[int]:
        """读取电荷与自旋多重度"""
        finds=re.search('(-?\d) (\d)',self.text)
        assert finds is not None,'未正确读取电荷与自旋多重度!'
        charge,multi=finds.groups()
        return int(charge),int(multi)