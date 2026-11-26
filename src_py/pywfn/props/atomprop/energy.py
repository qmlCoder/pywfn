"""
计算分子轨道内每个原子的能量
"""
from pywfn.base.mole import Mole

class Calculator:
    def __init__(self,mole:Mole) -> None:
        self.mole=mole