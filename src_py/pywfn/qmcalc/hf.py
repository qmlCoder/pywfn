from pywfn.base import Mole, Coefs
from pywfn import _core

import numpy as np


class RHF:
    def __init__(self, mole: Mole) -> None:
        self.core = _core.qmcalc.hf.RHF(mole.core)  # type: ignore

    def run(self, emat: np.ndarray):
        return self.core.run(emat)
