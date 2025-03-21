from pywfn.reader import Reader
from pywfn.data import consts

import numpy as np

class SdfReader(Reader):
    def __init__(self, path: str, cache: bool = False) -> None:
        super().__init__(path, cache)
        self.type='sdf'

    def get_atmSyms(self) -> list[str]:
        syms,xyzs=self.read_geom()
        return syms
    
    def get_atmXyzs(self) -> np.ndarray:
        syms,xyzs=self.read_geom()
        return xyzs

    def read_geom(self):
        finds=self.getline(3).split()
        natm=int(finds[0])
        lines=self.getlines(4,4+natm)
        xyzs=[]
        syms=[]
        for line in lines:
            x,y,z,s=line.split()[:4]
            xyzs.append([float(x),float(y),float(z)])
            syms.append(s)
        return syms,np.array(xyzs)/consts.Bohr

