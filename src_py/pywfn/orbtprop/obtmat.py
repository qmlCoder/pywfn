from pywfn import _core
from pywfn.base import Mole, Stm


class Calculator:
    def __init__(self, mole: Mole, atms: list[int] | None = None) -> None:
        self.mole = mole
        self.core = _core.orbtprop.obtmat.Calculator(mole.core, atms)  # type: ignore

    def pocv(
        self, dirs: dict[int, list[float]], keep_other_atm: bool, keep_other_sym: bool
    ):
        return self.core.pocv(dirs, keep_other_atm, keep_other_sym)

    def pi_pocv(self):
        return self.core.pi_pocv()

    def mocv(self, stms: dict[int, "Stm"], atos: list[int]):
        stms = {atm: stm.core for atm, stm in stms.items()}
        return self.core.mocv(stms, atos)

    def pi_mocv(self, dest: int = 0):
        return self.core.pi_mocv(dest)

    def ldao(self, stms: dict[int, "Stm"]):
        stms = {atm: stm.core for atm, stm in stms.items()}
        return self.core.ldao(stms)


def mocv(L, rcoefs, keeps):
    return _core.orbtprop.obtmat.mocv(L, rcoefs, keeps)  # type: ignore
