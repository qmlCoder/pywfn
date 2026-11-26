"""
得到分子的pi分子轨道
"""

from pywfn.base.mole import Mole
from pywfn.maths.mol import projCM
from pywfn.atomprop import direction
import numpy as np

class Calculator:
    def __init__(self,mole:Mole) -> None:
        self.mole=mole

    def project(self):
        dirCaler=direction.Calculator(self.mole)
        atms=[]
        dirs=[]
        for atom in self.mole.atoms:
            normal=dirCaler.normal_vector(atom.idx)
            if normal is None:continue
            atms.append(atom.idx)
            x,y,z=normal
            dirs.append([x,y,z])
        dirs=np.array(dirs)
        obts=self.mole.O_obts
        CMp=projCM(self.mole,obts,atms,dirs,False,False)
        return CMp