"""
计算分子的布居分析
"""

from pywfn import _core
from pywfn.base.mole import Mole


class Calculator:
    def __init__(self, mole: Mole) -> None:
        self.mole = mole
        self.core = _core.props.moleprop.population.Calculator(mole.core)  # type: ignore
