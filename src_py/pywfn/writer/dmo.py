"""
将分子对象保存为dmo文件
"""

from pywfn import _core
from pywfn.base.mole import Mole


class DmoWriter:
    def __init__(self) -> None:
        """dmo文件写出器"""
        self.core = _core.writer.DmoWriter()  # type: ignore

    @property
    def syms(self) -> list[str]:
        return self.core.syms

    @syms.setter
    def syms(self, syms: list[str]) -> None:
        self.core.syms = syms

    @property
    def xyzs(self) -> list[list[float]]:
        return self.core.xyzs

    @xyzs.setter
    def xyzs(self, xyzs: list[list[float]]) -> None:
        self.core.xyzs = xyzs

    @property
    def atos(self) -> list[tuple[int, int, str, list[float]]]:
        return self.core.atos

    @atos.setter
    def atos(self, atos: list[tuple[int, int, str, list[float]]]) -> None:
        self.core.atos = atos

    @property
    def BM(self) -> list[list[float]]:
        return self.core.BM

    @BM.setter
    def BM(self, bm: list[list[float]]) -> None:
        self.core.BM = bm

    @property
    def CM(self) -> list[list[float]]:
        return self.core.CM

    @CM.setter
    def CM(self, cm: list[list[float]]) -> None:
        self.core.CM = cm

    @property
    def es(self) -> list[float]:
        return self.core.es

    @es.setter
    def es(self, es: list[float]) -> None:
        self.core.es = es

    def save(self, path: str) -> None:
        """保存dmo文件"""
        self.core.save(path)

    def read_mole(self, mole: Mole) -> None:
        """从分子对象中读取数据"""
        self.core.read_mole(mole.core)
