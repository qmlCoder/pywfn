"""
定义mol文件的写出器
"""
from pywfn.base.mole import Mole
from pywfn.datas import temps
from pathlib import Path
import numpy as np
from pywfn import _core


class MolWriter:
    def __init__(self):
        self.core = _core.writer.MolWriter()  # type: ignore

    def read_mole(self, mole: Mole):
        self.core.read_mole(mole.core)

    def build(self, mole: Mole):
        return self.core.build(mole.core)

    def save(self, path: str):
        self.core.save(path)
