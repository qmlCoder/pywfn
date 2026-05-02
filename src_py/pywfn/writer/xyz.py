"""
将分子导出为xyz文件
"""

from pywfn.base.mole import Mole
from pywfn.datas import temps
from pywfn import config
import numpy as np
from pywfn import _core


class XyzWriter():  # type: ignore
    def __init__(self) -> None:
        self.core = _core.writer.XyzWriter()

    def read_mole(self, mole: Mole):
        self.core.read_mole(mole)

    def build(self) -> str:
        return self.core.build()

    def save(self, path: str):
        self.core.save(path)
