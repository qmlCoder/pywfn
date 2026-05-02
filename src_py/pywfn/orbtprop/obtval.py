from pywfn import _core
from pywfn.base import Mole


class Calculator:
    def __init__(self, mole: Mole) -> None:
        self.mole = mole
        self.core = _core.orbtprop.obtval.Calculator(mole.core)  # type: ignore

    def pi_eles_pocv(self) -> list[float]:
        return self.core.pi_eles_pocv()

    def pi_eles_ldao(self) -> list[float]:
        return self.core.pi_eles_ldao()
