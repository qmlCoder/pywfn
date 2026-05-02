from pywfn.base.mole import Mole
from pywfn.gridprop import dftgrid, density
from pywfn import _core

import numpy as np


class Calculator:
    def __init__(self, mole: Mole) -> None:
        self.mole = mole
