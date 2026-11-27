"""
将分子导出为xyz文件
"""

from pywfn import _core
from pywfn.base.mole import Mole


class XyzWriter:
    """标准 XYZ 坐标文件写出器"""

    def __init__(self) -> None:
        """xyz文件写出器"""
        self.core = _core.writer.XyzWriter()  # type: ignore

    @property
    def syms(self) -> list[str]:
        """原子符号列表"""
        return self.core.syms

    @syms.setter
    def syms(self, syms: list[str]) -> None:
        self.core.syms = syms

    @property
    def xyzs(self) -> list[list[float]]:
        """原子坐标列表"""
        return self.core.xyzs

    @xyzs.setter
    def xyzs(self, xyzs: list[list[float]]) -> None:
        self.core.xyzs = xyzs

    @property
    def title(self) -> str:
        """标题（文件的第二行）"""
        return self.core.title

    @title.setter
    def title(self, title: str) -> None:
        self.core.title = title

    def read_mole(self, mole: Mole) -> None:
        """从分子对象中读取数据

        Args:
            mole: 分子对象
        """
        self.core.read_mole(mole.core)

    def build(self) -> str:
        """生成xyz文件内容"""
        return self.core.build()

    def save(self, path: str) -> None:
        """保存xyz文件

        Args:
            path: 保存路径
        """
        self.core.save(path)
