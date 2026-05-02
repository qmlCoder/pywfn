from pywfn.base import Stm
from pywfn import _core
import numpy as np


class Calculator:
    def __init__(self, mole) -> None:
        self.mole = mole
        self.caler = _core.moleprop.orbital.Calculator(mole.core)  # type: ignore # 核心计算器



