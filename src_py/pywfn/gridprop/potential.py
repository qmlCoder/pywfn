from pywfn.base.mole import Mole
from pywfn.gridprop import density, dftgrid
import numpy as np
import json
from pathlib import Path
from pywfn import _core


class Calculator:
    def __init__(self, mole: Mole) -> None:
        self.mole = mole
        self.core = _core.gridprop.potential.Calculator(mole.core)  # type: ignore
