from pywfn.base import Mole
from pywfn import _core
from pywfn.base import Stm


class Calculator:
    def __init__(self, mole: Mole, atms: list[int] | None = None) -> None:
        self.mole = mole
        self.core = _core.orbtprop.obtmat.Calculator(
            mole.core, atms)  # type: ignore

    def pocv(self, dirs: dict[int, list[float]], keep_other_atm: bool, keep_other_sym: bool):
        return self.core.pocv(dirs, keep_other_atm, keep_other_sym)

    def pi_pocv(self):
        return self.core.pi_pocv()

    def mocv(self, mocvs: dict[int, "Mocv"], atos: list[int] | None = None):
        mocvs = {atm: mocv.core for atm, mocv in mocvs.items()}
        return self.core.mocv(mocvs, atos)

    def pi_mocv(self, mocvs: dict[int, "Mocv"]):
        return self.core.pi_mocv(mocvs)


class Mocv:
    def __init__(self, core: _core.orbtprop.obtmat.Mocv) -> None:  # type: ignore
        self.core = core  # type: ignore # 核心计算器

    @staticmethod
    def new(stm: Stm, skeep: list[int], pkeep: list[int], dkeep: list[int]):
        core = _core.orbtprop.obtmat.Mocv(stm.core, skeep, pkeep, dkeep)  # type: ignore
        return Mocv(core)
