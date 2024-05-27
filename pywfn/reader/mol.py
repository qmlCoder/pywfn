"""
mol文件读取器
"""
import numpy as np
from numpy import ndarray
from pywfn import reader
from functools import lru_cache
class MolReader(reader.Reader):
    def __init__(self, path: str) -> None:
        super().__init__(path)
        line3=self.getline(3)
        self.natm=int(line3[:3])
        print(self.natm)

    def get_coords(self) -> ndarray:
        return self.read_coords()[1]
    
    def get_symbols(self) -> list[str]:
        return self.read_coords()[0]
    
    @lru_cache
    def read_coords(self):
        """读取原子坐标"""
        xyzs=[]
        syms=[]
        for i in range(4,4+self.natm):
            line=self.getline(i)[:-1]
            x=line[ 0:10].strip()
            y=line[10:20].strip()
            z=line[20:30].strip()
            s=line[30:32].strip()
            xyzs.append([x,y,z])
            syms.append(s)
        return syms,np.array(xyzs,dtype=float)*1.889
