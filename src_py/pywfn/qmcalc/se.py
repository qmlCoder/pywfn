from pywfn.base import Mole, Coefs
from pywfn import _core


class RSE:
    def __init__(self, mole: Mole) -> None:
        self.core = _core.qmcalc.se.RSE(mole.core)  # type: ignore

    def run(self):
        return self.core.run()
