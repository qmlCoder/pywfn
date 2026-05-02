from pywfn import _core
from pywfn.base import Mole, Stm
from pywfn.atomprop import direction


class Int_mat:
    def __init__(self, mole: Mole) -> None:
        self.mole = mole
        self.core = _core.integ.Int_mat(mole.core)  # type: ignore

    def smat(self, form: str):
        return self.core.smat(form)

    def tmat(self, form: str):
        return self.core.tmat(form)

    def vmat(self, form: str):
        return self.core.vmat(form)

    def mat1e(self, form: str):
        return self.core.mat1e(form)

    def mat2e(self):
        return self.core.mat2e()

    def mat1e_lcs(self, stms: dict[int, Stm], cals: list[int]):
        stms = {atm: stm.core for atm, stm in stms.items()}
        return self.core.mat1e_lcs(stms, cals)
