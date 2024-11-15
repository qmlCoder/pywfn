"""
定义gjf的读取器，感觉可能也不会用得到
"""
from numpy import ndarray
from pywfn import reader
from pywfn.data.elements import elements
import re
import numpy as np
from functools import lru_cache


class GjfReader(reader.Reader):
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
        finds:list[str]=re.findall(r' ([A-Za-z\d]+) +(-?\d+.\d+) +(-?\d+.\d+) +(-?\d+.\d+)',self.text)
        symbols=[]
        coords=[]
        for s,x,y,z in finds:
            if s.isdigit():
                s=elements[int(s)].symbol
            symbols.append(s)
            coords.append([x,y,z])
        coords=np.array(coords,dtype=np.float32)*1.889
        return symbols,coords
    
    @lru_cache
    def read_multi(self)->tuple[int,int]:
        """读取电荷与自旋多重度"""
        finds=re.search(r'(-?\d) (\d)',self.text)
        assert finds is not None,'未正确读取电荷与自旋多重度!'
        charge,multi=finds.groups()
        return int(charge),int(multi)