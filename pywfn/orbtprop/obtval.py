from pywfn import core
from pywfn.base import Mole

class Calculator:
    def __init__(self,mole:Mole) -> None:
        self.mole=mole
        self.core=core.orbtprop.obtval.Calculator(mole.mole) # type: ignore

    def pi_ele(self)->list[float]:
        return self.core.pi_ele()