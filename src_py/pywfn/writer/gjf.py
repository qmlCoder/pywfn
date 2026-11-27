"""
将分子对象保存为gif文件
"""

from pywfn import _core
from pywfn.base.mole import Mole


class GjfWriter:
    """Gaussian 输入文件 (.gjf) 写出器"""

    def __init__(self):
        self.core = _core.writer.GjfWriter()  # type: ignore


    def set(self, syms: list[str], xyzs: list[list[float]], chk: str , job:str, charge: int, multip: int) -> None:
        """
        设置gjf文件的内容

        Args:
            syms: 原子符号列表
            xyzs: 原子坐标列表
            chk: checkpoint 文件名
            job: 任务类型，如 "opt"、"freq"
            charge: 分子总电荷
            multip: 自旋多重度
        """
        self.core.set(syms, xyzs, chk, job, charge, multip)

    @property
    def chk(self) -> str:
        """checkpoint 文件名"""
        return self.core.chk

    @chk.setter
    def chk(self, chk: str) -> None:
        self.core.chk = chk

    @property
    def job(self) -> str:
        """任务类型，如 "opt"、"freq" """
        return self.core.job

    @job.setter
    def job(self, job: str) -> None:
        self.core.job = job

    @property
    def charge(self) -> int:
        """分子总电荷"""
        return self.core.charge

    @charge.setter
    def charge(self, charge: int) -> None:
        self.core.charge = charge

    @property
    def multip(self) -> int:
        """自旋多重度"""
        return self.core.multip

    @multip.setter
    def multip(self, multip: int) -> None:
        self.core.multip = multip

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

    def build(self) -> str:
        """生成gjf文件内容"""
        return self.core.build()

    def save(self, path: str) -> None:
        """保存gjf文件"""
        return self.core.save(path)

    def read_mole(self, mole: Mole) -> None:
        """从分子对象中读取数据"""
        self.core.read_mole(mole.core)
