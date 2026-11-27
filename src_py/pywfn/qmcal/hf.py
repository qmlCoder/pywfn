from pywfn import _core
from pywfn.base import Mole


class RHF:
    def __init__(self, mole: Mole) -> None:
        self.core = _core.qmcal.hf.RHF(mole.core)  # type: ignore

    def run(self) -> None:
        """运行Hartree-Fock计算"""
        self.core.run()
