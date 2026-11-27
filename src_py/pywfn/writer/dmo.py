"""
将分子对象保存为dmo文件
"""

from pywfn import _core
from pywfn.base.mole import Mole


class DmoWriter:
    """DMol3 文件写出器"""

    def __init__(self) -> None:
        """dmo文件写出器"""
        self.core = _core.writer.DmoWriter()  # type: ignore

    @staticmethod
    def blank() -> "DmoWriter":
        """创建一个空的写出器，与 `DmoWriter()` 等价

        Returns:
            DmoWriter: 空的 dmo 写出器
        """
        writer = DmoWriter.__new__(DmoWriter)
        writer.core = _core.writer.DmoWriter.blank()  # type: ignore
        return writer

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
    def atos(self) -> list[tuple[int, int, str, list[float]]]:
        """原子轨道列表 (原子索引, 壳层索引, 符号, 坐标)"""
        return self.core.atos

    @atos.setter
    def atos(self, atos: list[tuple[int, int, str, list[float]]]) -> None:
        self.core.atos = atos

    @property
    def BM(self) -> list[list[float]]:
        """键连矩阵"""
        return self.core.BM

    @BM.setter
    def BM(self, bm: list[list[float]]) -> None:
        self.core.BM = bm

    @property
    def CM(self) -> list[list[float]]:
        """轨道系数矩阵"""
        return self.core.CM

    @CM.setter
    def CM(self, cm: list[list[float]]) -> None:
        self.core.CM = cm

    @property
    def es(self) -> list[float]:
        """轨道能量列表"""
        return self.core.es

    @es.setter
    def es(self, es: list[float]) -> None:
        self.core.es = es

    def save(self, path: str) -> None:
        """保存dmo文件

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
