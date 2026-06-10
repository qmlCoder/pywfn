"""
将分子对象保存为gif文件
"""

from pywfn import _core
from pywfn.base.mole import Mole


class GjfWriter:
    def __init__(self):
        self.core = _core.writer.GjfWriter()  # type: ignore

    @property
    def syms(self) -> list[str]:
        return self.core.syms

    @syms.setter
    def syms(self, syms: list[str]):
        self.core.syms = syms

    @property
    def xyzs(self) -> list[list[float]]:
        return self.core.xyzs

    @xyzs.setter
    def xyzs(self, xyzs: list[list[float]]):
        self.core.xyzs = xyzs

    @property
    def chk(self) -> str:
        return self.core.chk

    @chk.setter
    def chk(self, chk: str):
        self.core.chk = chk

    @property
    def job(self) -> str:
        return self.core.job

    @job.setter
    def job(self, job: str):
        self.core.job = job

    @property
    def charge(self) -> int:
        return self.core.charge

    @charge.setter
    def charge(self, charge: int):
        self.core.charge = charge

    @property
    def multip(self) -> int:
        return self.core.multip

    @multip.setter
    def multip(self, multip: int):
        self.core.multip = multip

    def from_mole(self, mole: Mole):
        self.core.from_mole(mole.core)

    def build(self) -> str:
        return self.core.build()

    def save(self, path: str):
        return self.core.save(path)
