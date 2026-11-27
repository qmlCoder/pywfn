"""
各种波函数文件的读取器

本模块对应 rswfn 的 `reader` 模块（Rust 侧由 `_core.reader` 提供）。
每种文件格式放在单独的子模块里（如 `pywfn.reader.mfn`），
各格式共用的数据类（`FileType`/`BasData`/`Marks`/`Vibe`/`GauJob`/`Nics`）
与工具函数（`read_lines`/`charge_number`/`get_fileType`/`parse_job`/
`shl_basis`/`shl_forms`）放在本文件里，最后统一从这里导出。
"""

from enum import IntEnum
from pathlib import Path
from typing import overload

from pywfn import _core
from pywfn.base import Atoms, Basis, Mole, Stm


class Nics:
    """NICS（核独立化学位移）数据"""

    def __init__(self, core: "_core.reader.Nics") -> None:  # type: ignore
        """包装 `_core` 的 Nics 对象

        Args:
            core: `_core.reader.Nics` 对象
        """
        self.core = core

    @property
    def iso(self) -> float:
        """各向同性化学位移"""
        return self.core.iso

    @property
    def ani(self) -> float:
        """各向异性化学位移"""
        return self.core.ani

    @property
    def mat(self) -> list[float]:
        """屏蔽张量（长度 9 的列表）"""
        return self.core.mat

    @property
    def eig(self) -> list[float]:
        """屏蔽张量的三个本征值"""
        return self.core.eig


class FileType(IntEnum):
    """支持读取的文件类型

    与 rswfn 的 `reader::FileType` 枚举、以及 `_core.reader.FileType` 一一对应。
    """

    Gjf = 0  # Gaussian 输入文件 (.gjf)
    Log = 1  # Gaussian 输出文件 (.log / .out)
    Fch = 2  # Gaussian 格式化检查点文件 (.fch / .fchk)
    Mol = 3  # MDL 分子文件 (.mol)
    Sdf = 4  # MDL 结构数据文件 (.sdf)
    Xyz = 5  # 标准 XYZ 坐标文件 (.xyz)
    Mde = 6  # Molden 波函数文件 (.molden)
    Mfn = 7  # Multiwfn 波函数文件 (.mwfn)
    Dmo = 8  # DMol3 输出文件 (.dmo / .dmi)
    Inp = 9  # 通用输入文件 (.inp)
    Cub = 10  # 立方体文件 (.cub)


# Python 枚举与 `_core` 枚举之间的对应关系，只在模块加载时建立一次
_CORE_FILETYPE: dict[FileType, "_core.reader.FileType"] = {
    FileType.Gjf: _core.reader.FileType.Gjf,  # type: ignore
    FileType.Log: _core.reader.FileType.Log,  # type: ignore
    FileType.Fch: _core.reader.FileType.Fch,  # type: ignore
    FileType.Mol: _core.reader.FileType.Mol,  # type: ignore
    FileType.Sdf: _core.reader.FileType.Sdf,  # type: ignore
    FileType.Xyz: _core.reader.FileType.Xyz,  # type: ignore
    FileType.Mde: _core.reader.FileType.Mde,  # type: ignore
    FileType.Mfn: _core.reader.FileType.Mfn,  # type: ignore
    FileType.Dmo: _core.reader.FileType.Dmo,  # type: ignore
    FileType.Inp: _core.reader.FileType.Inp,  # type: ignore
    FileType.Cub: _core.reader.FileType.Cub,  # type: ignore
}


def _file_type_from_core(core_type: "_core.reader.FileType") -> FileType:
    """把 `_core` 的文件类型枚举转换成 Python 枚举

    Args:
        core_type: `_core.reader.FileType` 枚举值

    Returns:
        FileType: 对应的 Python 枚举值

    Raises:
        ValueError: `_core` 中出现了 Python 层未收录的枚举值
    """
    for ftype, ctype in _CORE_FILETYPE.items():
        if ctype == core_type:
            return ftype
    raise ValueError(f"未知的文件类型枚举值: {core_type}")


class BasData:
    """基函数的原始数据

    对应 rswfn 中的类型别名 `BasData = (usize, usize, usize, f64, f64)`，
    即 (原子索引, 壳层索引, 角动量, 高斯指数, 收缩系数)。
    """

    def __init__(self, atm: int, shl: int, ang: int, alp: float, coe: float) -> None:
        """创建一个基函数原始数据对象

        Args:
            atm: 所属原子索引
            shl: 所属壳层索引
            ang: 角动量量子数
            alp: 高斯指数
            coe: 收缩系数
        """
        self.core = _core.reader.BasData(atm, shl, ang, alp, coe)  # type: ignore

    @staticmethod
    def from_core(core: "_core.reader.BasData") -> "BasData":
        """把 `_core` 的 BasData 包装成 Python 对象

        Args:
            core: `_core.reader.BasData` 对象

        Returns:
            BasData: 包装后的对象
        """
        data = BasData.__new__(BasData)
        data.core = core
        return data

    @property
    def atm(self) -> int:
        """所属原子索引"""
        return self.core.atm

    @property
    def shl(self) -> int:
        """所属壳层索引"""
        return self.core.shl

    @property
    def ang(self) -> int:
        """角动量量子数"""
        return self.core.ang

    @property
    def alp(self) -> float:
        """高斯指数"""
        return self.core.alp

    @property
    def coe(self) -> float:
        """收缩系数"""
        return self.core.coe

    def as_tuple(self) -> tuple[int, int, int, float, float]:
        """转换为元组

        Returns:
            tuple[int, int, int, float, float]: (atm, shl, ang, alp, coe)
        """
        return self.core.as_tuple()

    def __repr__(self) -> str:
        return self.core.__repr__()


class Marks:
    """文件中的关键字标记表

    记录每个关键字（如 "NAtoms"）在文件行列表中的行号，用于快速定位数据块。
    """

    def __init__(self) -> None:
        """创建一个空的标记表"""
        self.core = _core.reader.Marks()  # type: ignore

    @property
    def marks(self) -> dict[str, int]:
        """全部标记及其行号"""
        return self.core.marks

    def add(self, marker: str, number: int) -> None:
        """添加一个标记

        Args:
            marker: 标记名称
            number: 标记所在的行号
        """
        self.core.add(marker, number)


class Vibe:
    """一个正则振动模式

    包含振动频率、约化质量、力常数、红外强度以及各原子的振动位移。
    """

    def __init__(self, core: "_core.reader.Vibe") -> None:  # type: ignore
        """包装 `_core` 的 Vibe 对象

        Args:
            core: `_core.reader.Vibe` 对象
        """
        self.core = core

    @property
    def freq(self) -> float:
        """振动频率"""
        return self.core.freq

    @property
    def redm(self) -> float:
        """约化质量"""
        return self.core.redm

    @property
    def frcc(self) -> float:
        """力常数"""
        return self.core.frcc

    @property
    def irin(self) -> float:
        """红外强度"""
        return self.core.irin

    @property
    def vibs(self) -> list[list[float]]:
        """每个原子的振动位移向量"""
        return self.core.vibs

    def get_max(self) -> tuple[int, list[float], float]:
        """获取振动幅度最大的原子及其位移方向

        Returns:
            tuple[int, list[float], float]: (原子索引, 该原子的位移向量, 该原子的位移大小)
        """
        return self.core.get_max()


class GauJob:
    """Gaussian 任务信息（从 # 行解析得到）"""

    def __init__(self, core: "_core.reader.GauJob") -> None:  # type: ignore
        """包装 `_core` 的 GauJob 对象

        Args:
            core: `_core.reader.GauJob` 对象
        """
        self.core = core

    @property
    def job_type(self) -> str:
        """任务类型（如 opt、freq、sp 等）"""
        return self.core.job_type

    @property
    def functional(self) -> str | None:
        """泛函名称（如 B3LYP），未指定时为 None"""
        return self.core.functional

    @property
    def basis_set(self) -> str | None:
        """基组名称（如 6-31G(d)），未指定时为 None"""
        return self.core.basis_set

    def __repr__(self) -> str:
        return (
            f"GauJob(job_type={self.job_type!r}, functional={self.functional!r}, "
            f"basis_set={self.basis_set!r})"
        )


def _wrap_atoms(core_atoms: "_core.base.Atoms") -> Atoms:
    """把 `_core` 的 Atoms 包装成 Python 的 Atoms

    由于 Python 的 `Atoms` 需要绑定一个分子对象，这里用同一份原子信息
    构建一个只含原子的 `Mole`（不含基组与轨道系数）。

    Args:
        core_atoms: `_core.base.Atoms` 对象

    Returns:
        Atoms: 包装后的原子集合
    """
    core_mole = _core.base.Mole.build(core_atoms, None, None)  # type: ignore
    return Atoms(core_atoms, Mole(core_mole))


def _wrap_stm(core_stm: "_core.base.Stm") -> Stm:
    """把 `_core` 的 Stm 包装成 Python 的 Stm

    Args:
        core_stm: `_core.base.Stm` 对象

    Returns:
        Stm: 包装后的局域坐标系
    """
    return Stm(core_stm)


def _wrap_bas_data(core_data: "_core.reader.BasData") -> BasData:
    """把 `_core` 的 BasData 包装成 Python 的 BasData

    Args:
        core_data: `_core.reader.BasData` 对象

    Returns:
        BasData: 包装后的基函数原始数据
    """
    return BasData.from_core(core_data)


def _wrap_vibe(core_vibe: "_core.reader.Vibe") -> Vibe:
    """把 `_core` 的 Vibe 包装成 Python 的 Vibe

    Args:
        core_vibe: `_core.reader.Vibe` 对象

    Returns:
        Vibe: 包装后的振动模式
    """
    return Vibe(core_vibe)


def _wrap_gau_job(core_job: "_core.reader.GauJob") -> GauJob:
    """把 `_core` 的 GauJob 包装成 Python 的 GauJob

    Args:
        core_job: `_core.reader.GauJob` 对象

    Returns:
        GauJob: 包装后的任务信息
    """
    return GauJob(core_job)


def read_lines(path: str) -> list[str]:
    """读取文件的全部行

    Args:
        path: 文件路径

    Returns:
        list[str]: 文件每一行组成的字符串列表

    Raises:
        FileNotFoundError: 文件不存在
    """
    if not Path(path).exists():
        raise FileNotFoundError(f"文件不存在: {path}")
    return _core.reader.read_lines(path)  # type: ignore


def charge_number(charge: int, multip: int, total: int) -> tuple[int, int]:
    """由体系电荷、自旋多重度和总电子数计算 alpha/beta 电子数

    Args:
        charge: 体系电荷
        multip: 自旋多重度
        total: 总电子数

    Returns:
        tuple[int, int]: (alpha 电子数, beta 电子数)
    """
    nela, nelb = _core.reader.charge_number(charge, multip, total)  # type: ignore
    return nela, nelb


def get_fileType(path: str) -> FileType:
    """根据文件后缀判断文件类型

    Args:
        path: 文件路径

    Returns:
        FileType: 对应的文件类型枚举

    Raises:
        ValueError: 无法识别的文件类型
    """
    try:
        core_type = _core.reader.get_fileType(path)  # type: ignore
    except BaseException:  # _core 遇到未知后缀会 panic
        raise ValueError(f"无法识别的文件类型: {path}") from None
    return _file_type_from_core(core_type)


def parse_job(line: str) -> GauJob | None:
    """解析高斯作业关键词行（形如 `#p b3lyp/6-31g(d) opt freq`）

    Args:
        line: Gaussian 输入文件中的路由行（以 # 开头的那一行）

    Returns:
        GauJob | None: 解析成功时返回任务信息，否则返回 None
    """
    core_job = _core.reader.parse_job(line)  # type: ignore
    return None if core_job is None else _wrap_gau_job(core_job)


# 内置基组库：{基组名: 基组定义原文}，只在模块加载时建立一次
_BASIS_SETS: dict[str, str] = {
    name: text for name, text in _core.rsbse.BASIS_SETS  # type: ignore
}


def _split_basis_text(text: str) -> dict[str, str]:
    """把基组定义原文按元素拆开

    基组原文是 Gaussian 格式，每个元素以 `C     0` 这样的一行开头。

    Args:
        text: 基组定义原文

    Returns:
        dict[str, str]: {元素符号: 该元素的基组定义文本}
    """
    blocks: dict[str, str] = {}
    symbol: str | None = None
    lines: list[str] = []
    for line in text.splitlines():
        parts = line.split()
        is_head = len(parts) == 2 and parts[1] == "0" and parts[0].isalpha()
        if is_head:
            if symbol is not None:
                blocks[symbol] = "\n".join(lines).strip()
            symbol = parts[0]
            lines = [line]
        elif symbol is not None:
            lines.append(line)
    if symbol is not None:
        blocks[symbol] = "\n".join(lines).strip()
    return blocks


def shl_basis(datas: list[BasData], forms: list[str]) -> Basis:
    """由基函数原始数据和壳层形式构建基组

    Args:
        datas: 每个基函数的原始数据
        forms: 每个壳层的基函数形式（"sph" 或 "car"）

    Returns:
        Basis: 构建好的基组对象
    """
    core_datas = [data.core for data in datas]
    core_basis = _core.reader.shl_basis(core_datas, forms)  # type: ignore
    # 基组本身不携带分子，这里绑定一个只有基组的空分子，便于占位
    empty_atoms = _core.base.Atoms([], [], [])  # type: ignore
    core_mole = _core.base.Mole.build(empty_atoms, core_basis, None)  # type: ignore
    return Basis(core_basis, Mole(core_mole))


def shl_forms(
    ato_atms: list[int], ato_shls: list[int], ato_syms: list[str]
) -> list[str]:
    """确定每个壳层的基函数形式（球谐或笛卡尔）

    Args:
        ato_atms: 每个基函数所属的原子索引
        ato_shls: 每个基函数所属的壳层索引
        ato_syms: 每个基函数的符号

    Returns:
        list[str]: 每个壳层对应的形式（"sph" 或 "car"）
    """
    return _core.reader.shl_forms(ato_atms, ato_shls, ato_syms)  # type: ignore


from pywfn.reader.cub import CubReader  # noqa: E402
from pywfn.reader.dmo import DmoReader  # noqa: E402
from pywfn.reader.fch import FchReader  # noqa: E402
from pywfn.reader.gjf import GjfReader  # noqa: E402
from pywfn.reader.inp import InpReader  # noqa: E402
from pywfn.reader.log import LogReader  # noqa: E402
from pywfn.reader.mde import MdeReader  # noqa: E402
from pywfn.reader.mfn import MfnReader  # noqa: E402
from pywfn.reader.mol import MolReader  # noqa: E402
from pywfn.reader.non import NonReader  # noqa: E402
from pywfn.reader.out import OutReader  # noqa: E402
from pywfn.reader.sdf import SdfReader  # noqa: E402
from pywfn.reader.xyz import XyzReader  # noqa: E402

__all__ = [
    "BasData",
    "CubReader",
    "DmoReader",
    "FchReader",
    "FileType",
    "GauJob",
    "GjfReader",
    "InpReader",
    "LogReader",
    "Marks",
    "MdeReader",
    "MfnReader",
    "MolReader",
    "Nics",
    "NonReader",
    "OutReader",
    "SdfReader",
    "Vibe",
    "XyzReader",
    "charge_number",
    "get_fileType",
    "get_reader",
    "parse_job",
    "read_lines",
    "shl_basis",
    "shl_forms",
    "supports",
]

# 与 rswfn 的 `reader::get_fileType` 保持一致
# （`.mfn` 是 Multiwfn 波函数文件的常用后缀，rswfn 的 get_fileType 只认 `.mwfn`，
#   这里在 get_reader 里额外补上，方便直接读取 .mfn 文件）
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
    ".mwfn",
    ".mfn",
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
    ".mwfn": MfnReader,
    ".mfn": MfnReader,
    ".dmi": DmoReader,
    ".dmo": DmoReader,
    ".inp": InpReader,
    ".cub": CubReader,
}


def get_reader(
    path: str,
) -> "FchReader | LogReader | GjfReader | MolReader | SdfReader | XyzReader | MdeReader | MfnReader | DmoReader | InpReader | CubReader":
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
