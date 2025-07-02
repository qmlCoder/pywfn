"""
得到分子的pi分子轨道
"""

from pywfn.base import Mole
from pywfn.maths.mol import projCM
from pywfn.atomprop import direction
import numpy as np

class Calculator:
    def __init__(self,mol:Mole) -> None:
        self.mol=mol

    def project(self):
        dirCaler=direction.Calculator(self.mol)
        atms=[]
        dirs=[]
        for atom in self.mol.atoms:
            normal=dirCaler.normal(atom.idx)
            if normal is None:continue
            atms.append(atom.idx)
            x,y,z=normal
            dirs.append([x,y,z])
        dirs=np.array(dirs)
        obts=self.mol.O_obts
        CMp=projCM(self.mol,obts,atms,dirs,False,False)
        return CMp