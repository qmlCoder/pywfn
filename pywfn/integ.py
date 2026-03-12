from pywfn import core
from pywfn.base import Mole,Stm
from pywfn.atomprop import direction

class Int_mat:
    def __init__(self,mole:Mole) -> None:
        self.mole=mole
        self.core=core.integ.Int_mat(mole.mole) # type: ignore

    def smat(self):
        return self.core.smat()
    
    def tmat(self):
        return self.core.tmat()
    
    def vmat(self):
        return self.core.vmat()
    
    def mat1e(self):
        return self.core.mat1e()
    
    def mat2e(self):
        return self.core.mat2e()
    
    def mat1e_lcs(self,stms:list[Stm]|None=None,level:int=2):
        natm=self.mole.syms.__len__()
        if stms is None:
            dir_caler=direction.Calculator(self.mole)
            stms=[dir_caler.LCS(atm) for atm in range(natm)]
        stms=[stm.core for stm in stms]
        return self.core.mat1e_lcs(stms,level)
