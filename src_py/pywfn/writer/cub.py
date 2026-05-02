"""
用来生成格点数据并保存cube文件
cube文件使用的是玻尔半径B
分子的坐标使用的是埃米A
A=B/1.889
B=A*1.889
计算的时候使用A，写文件的时候将坐标和格点转为B
"""
from numpy._core.numeric import ndarray
from numpy.typing import NDArray
from pywfn.base.mole import Mole
from pywfn import _core


class CubWriter:
    def __init__(self) -> None:
        """cube文件导出器"""
        self.core = _core.writer.CubWriter()  # type: ignore

    def set(self,
            title: str | None = None,  # 标题
            syms: list[str] | None = None,  # 元素符号
            xyzs: NDArray[float] | None = None,  # 原子坐标
            obts: list[int] | None = None,  # 轨道序数
            pos0: list[float] | None = None,  # 起始坐标
            size: list[int] | None = None,  # 格点数量
            step: list[float] | None = None,  # 格点步长
            vals: NDArray[float] | None = None  # 格点数值[ngrid,nobt]
            ):
        assert len(obts)==vals.shape[0], "obts和vals的维度不匹配"
        self.core.set(title, syms, xyzs, obts, pos0, size, step, vals)

    def read_mole(self, mole: Mole):
        """从分子对象中读取数据"""
        self.core.read_mole(mole.core)

    def save(self,path:str):
        self.core.save(path)
