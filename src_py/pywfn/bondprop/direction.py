"""
计算与键相关的方向
"""
from pywfn.base.mole import Mole
from pywfn import _core
import numpy as np


class Calculator:
    def __init__(self,mole:Mole) -> None:
        self.core=_core.bondprop.direction.Calculator(mole.core) # type: ignore

    def verts(self,ai:int,aj:int):
        """垂直于某根键的方向"""
        dirs=self.core.verts(ai,aj)
        return np.array(dirs)
        