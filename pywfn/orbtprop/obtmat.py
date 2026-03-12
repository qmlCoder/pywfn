from pywfn.base import Mole,Stm
from pywfn import core
from pywfn.atomprop import direction

class Calculator:
    def __init__(self,mole:Mole) -> None:
        self.mole=mole
        self.core=core.orbtprop.obtmat.Calculator(mole.mole) # type: ignore

    
    def ldao(self,stms:list[Stm]|None=None):
        """局部坐标系的分子轨道系数矩阵"""
        natm=self.mole.syms.__len__()
        if stms is None:
            dir_caler=direction.Calculator(self.mole)
            stms=[dir_caler.LCS(atm) for atm in range(natm)]
        stms=[stm.core for stm in stms]
        return self.core.ldao(stms)
