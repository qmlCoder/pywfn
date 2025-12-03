from pywfn.base import Mole
from pywfn import core

class Calculator:
    def __init__(self, mole:Mole):
        self.mole = mole
        self.caler=core.moleprop.newmole.Calculator(mole.mole) # type: ignore

    def mole_pi(self)->Mole:
        mole_core=self.caler.mole_pi()
        newmole=self.mole.clone()
        newmole.mole=mole_core
        return newmole