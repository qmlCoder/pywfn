from pywfn.base import Mole
from pywfn import _core


class Calculator:
    def __init__(self, mole: Mole):
        self.mole = mole
        self.caler = _core.moleprop.newmole.Calculator(mole.core)  # type: ignore
