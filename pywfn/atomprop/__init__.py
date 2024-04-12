"""
该模块用来计算原子属性
Mulliken电荷
"""
from pywfn.utils import printer

import numpy as np

class AtomCaler:
    def __init__(self) -> None:
        self.logTip:str=''

    def resStr(self)->np.ndarray:pass

    def logRes(self):
        printer.info(self.logTip)
        printer.res(self.resStr())