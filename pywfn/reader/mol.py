"""
mol文件读取器
"""
import numpy as np
from numpy import ndarray
from pywfn import reader
from pywfn.base.geome import Geome
from functools import lru_cache


class MolReader(reader.Reader):
    def __init__(self, path: str) -> None:
        super().__init__(path)
        self.type='mol'
        line3=self.getline(3)
        self.natm=int(line3[:3])
        print(self.natm)

    def get_geome(self) -> Geome:
        syms,xyzs=self.read_geome()
        return Geome().build(syms,xyzs)
    
    @lru_cache
    def read_geome(self):
        """读取原子坐标"""
        xyzs=[]
        syms:list[str]=[]
        for i in range(4,4+self.natm):
            line=self.getline(i)[:-1]
            x=line[ 0:10].strip()
            y=line[10:20].strip()
            z=line[20:30].strip()
            s=line[30:32].strip()
            xyzs.append([x,y,z])
            syms.append(s)
        return syms,np.array(xyzs,dtype=float)*1.889
