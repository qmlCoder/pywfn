"""
将分子对象保存为fch文件
"""

from pywfn import _core
from pywfn.base.mole import Mole


class FchWriter:
    def __init__(self) -> None:
        """fch文件写出器"""
        self.core = _core.writer.FchWriter()  # type: ignore

    @staticmethod
    def from_mole(mole: Mole) -> "FchWriter":
        """从分子对象构建写出器"""
        writer = FchWriter.__new__(FchWriter)
        writer.core = _core.writer.FchWriter.from_mole(mole.core)  # type: ignore
        return writer

    @property
    def charge(self) -> int:
        return self.core.charge

    @charge.setter
    def charge(self, charge: int) -> None:
        self.core.charge = charge

    @property
    def multip(self) -> int:
        return self.core.multip

    @multip.setter
    def multip(self, multip: int) -> None:
        self.core.multip = multip

    def build(self) -> str:
        """生成fch文件内容"""
        return self.core.build()

    def save(self, path: str) -> None:
        """保存fch文件"""
        self.core.save(path)

    def read_mole(self, mole: Mole) -> None:
        """从分子对象中读取数据"""
        self.core.read_mole(mole.core)
