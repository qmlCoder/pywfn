"""
将分子对象保存为gif文件
"""

from pywfn import _core
from pywfn.base.mole import Mole


class GjfWriter:
    def __init__(self):
        self.core = _core.writer.GjfWriter()  # type: ignore


    def set(self, syms: list[str], xyzs: list[list[float]], chk: str , job:str, charge: int, multip: int):
        """
        设置gjf文件的内容
        """
        self.core.set(syms, xyzs, chk, job, charge, multip)

    def build(self) -> str:
        """生成gjf文件内容"""
        return self.core.build()

    def save(self, path: str) -> None:
        """保存gjf文件"""
        return self.core.save(path)

    def read_mole(self, mole: Mole) -> None:
        """从分子对象中读取数据"""
        self.core.read_mole(mole.core)
