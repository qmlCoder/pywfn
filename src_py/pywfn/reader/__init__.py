class Nics:
    def __init__(self, core):
        self.core = core

    @property
    def iso(self) -> float:
        return self.core.iso

    @property
    def ani(self) -> float:
        return self.core.ani

    @property
    def mat(self) -> list[float]:
        return self.core.mat

    @property
    def eig(self) -> list[float]:
        return self.core.eig

from pathlib import Path

from pywfn.reader.cub import CubReader
from pywfn.reader.dmo import DmoReader
from pywfn.reader.fch import FchReader
from pywfn.reader.gjf import GjfReader
from pywfn.reader.inp import InpReader
from pywfn.reader.log import LogReader
from pywfn.reader.mde import MdeReader
from pywfn.reader.mol import MolReader
from pywfn.reader.non import NonReader
from pywfn.reader.out import OutReader
from pywfn.reader.sdf import SdfReader
from pywfn.reader.xyz import XyzReader

__all__ = [
    "CubReader",
    "DmoReader",
    "FchReader",
    "GjfReader",
    "InpReader",
    "LogReader",
    "MdeReader",
    "MolReader",
    "NonReader",
    "OutReader",
    "SdfReader",
    "XyzReader",
    "get_reader",
    "supports",
]

# from pywfn.reader.wfn import WfnReader

# 与 rswfn 的 `reader::get_fileType` 保持一致
supports = [
    ".gjf",
    ".log",
    ".out",
    ".fch",
    ".fchk",
    ".mol",
    ".sdf",
    ".xyz",
    ".molden",
    ".dmi",
    ".dmo",
    ".inp",
    ".cub",
]  # 支持的文件类型

# 后缀到读取器的映射，后缀规则与 rswfn 的 `get_fileType` 一致
_READERS = {
    ".gjf": GjfReader,
    ".log": LogReader,
    ".out": LogReader,
    ".fch": FchReader,
    ".fchk": FchReader,
    ".mol": MolReader,
    ".sdf": SdfReader,
    ".xyz": XyzReader,
    ".molden": MdeReader,
    ".dmi": DmoReader,
    ".dmo": DmoReader,
    ".inp": InpReader,
    ".cub": CubReader,
}


def get_reader(
    path: str,
) -> "FchReader | LogReader | GjfReader | MolReader | SdfReader | XyzReader | MdeReader | DmoReader | InpReader | CubReader":
    """根据输入文件的类型自动判断应该使用哪个读取器

    Args:
        path: 文件路径

    Returns:
        对应文件类型的读取器实例

    Raises:
        ValueError: 不支持的文件类型
    """
    suffix = Path(path).suffix.lower()
    if suffix not in _READERS:
        raise ValueError(f"不支持的文件类型: {suffix}，已支持: {supports}")
    return _READERS[suffix].from_path(path)  # type: ignore
