"""
将分子对象保存为 mfn（Multiwfn 波函数）文件
"""

from pywfn import _core
from pywfn.base.mole import Mole


class MfnWriter:
    """Multiwfn 波函数文件 (.mwfn) 写出器"""

    def __init__(self) -> None:
        """mfn文件写出器"""
        self.core = _core.writer.MfnWriter()  # type: ignore

    @staticmethod
    def from_mole(mole: Mole) -> "MfnWriter":
        """从分子对象构建写出器

        Args:
            mole: 分子对象

        Returns:
            MfnWriter: 已填好分子信息的写出器
        """
        writer = MfnWriter.__new__(MfnWriter)
        writer.core = _core.writer.MfnWriter.from_mole(mole.core)  # type: ignore
        return writer

    @property
    def mole(self) -> Mole:
        """分子对象"""
        return Mole(self.core.mole)

    @mole.setter
    def mole(self, mole: Mole) -> None:
        self.core.mole = mole.core

    @property
    def energy(self) -> float:
        """分子能量"""
        return self.core.energy

    @energy.setter
    def energy(self, energy: float) -> None:
        self.core.energy = energy

    def build(self) -> str:
        """生成mfn文件内容

        Returns:
            str: 格式化的 mfn 文件字符串
        """
        return self.core.build()

    def save(self, path: str) -> None:
        """保存mfn文件

        Args:
            path: 保存路径
        """
        self.core.save(path)

    def read_mole(self, mole: Mole) -> None:
        """从分子对象中读取数据

        Args:
            mole: 分子对象
        """
        self.core.read_mole(mole.core)
