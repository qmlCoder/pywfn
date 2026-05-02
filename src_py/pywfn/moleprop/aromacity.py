"""
计算分子的芳香性
使用pi键级的标准差表示
"""
from pywfn.base.mole import Mole
from pywfn.atomprop import charge
from pywfn.bondprop import order as orderProp
import numpy as np
from pywfn import _core
from pywfn.datas import consts


class Calculator:
    def __init__(self, mole: Mole):
        self.mole = mole
        self.ratio = 0.5
        self.caler = _core.moleprop.aromacity.Calculator(mole.core)  # type: ignore

    def PISD(self, rings: list[list[int]] | None = None):  # 版本1 直接用键级标准差
        return self.caler.PISD(rings)
