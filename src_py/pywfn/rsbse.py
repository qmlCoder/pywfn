"""
基组读取（Basis Set Exchange）

对应 rswfn 的 `rsbse` 模块，负责从基组文件、基组文本或 Basis Set Exchange
网站读取基组，并转换为 `pywfn.base.Basis` 对象。

内置的基组别名表见 `BASIS_SETS`，可用的基组名见 `supported_names`。
"""

from pywfn import _core
from pywfn.base import Atoms, Basis

BASIS_SETS: list[tuple[str, str]] = _core.rsbse.BASIS_SETS  # type: ignore
"""内置基组别名表，每个元素为 (基组名, 基组文本)"""


def from_file(path: str, atoms: Atoms) -> Basis | None:
    """从基组文件读取基组

    把 `path` 指定的基组文件内容解析为 Basis，只保留 `atoms` 中出现的元素。

    Args:
        path: 基组文件路径（如从 BSE 下载的 .gbs 文件）
        atoms: 原子集合对象，用于确定需要哪些元素的基组

    Returns:
        Basis | None: 解析失败时返回 None
    """
    core = _core.rsbse.from_file(path, atoms.core)  # type: ignore
    if core is None:
        return None
    return Basis(core, atoms.mole)


def from_name(name: str, atoms: Atoms) -> Basis | None:
    """按名称从 Basis Set Exchange 下载并解析基组

    名称先通过 `BASIS_SETS` 中的别名表匹配，找不到时按 BSE 上的名称下载。
    该函数会访问网络且不设置超时，可能耗时较长甚至失败。

    Args:
        name: 基组名称，如 "def2-SVP"
        atoms: 原子集合对象，用于确定需要哪些元素的基组

    Returns:
        Basis | None: 名称不存在或下载失败时返回 None
    """
    core = _core.rsbse.from_name(name, atoms.core)  # type: ignore
    if core is None:
        return None
    return Basis(core, atoms.mole)


def from_text(text: str, atoms: Atoms) -> Basis | None:
    """从基组文本内容解析基组

    Args:
        text: 基组文本内容（Gaussian 格式）
        atoms: 原子集合对象，用于确定需要哪些元素的基组

    Returns:
        Basis | None: 解析失败时返回 None
    """
    core = _core.rsbse.from_text(text, atoms.core)  # type: ignore
    if core is None:
        return None
    return Basis(core, atoms.mole)


def read_basis(lines: list[str]) -> list[tuple[int, int, int, float, float]] | None:
    """从基组文本的逐行内容中解析壳层数据

    Args:
        lines: 基组文本按行拆分后的列表

    Returns:
        list[tuple[int, int, int, float, float]] | None: 解析成功时返回列表，
        每个元素为 (原子序号, 壳层序号, 角动量, 高斯指数, 收缩系数)；解析失败时返回 None
    """
    return _core.rsbse.read_basis(lines)  # type: ignore


def supported_names() -> list[str]:
    """返回所有支持的基组名称

    Returns:
        list[str]: 基组名称列表，可直接传给 `from_name`
    """
    return _core.rsbse.supported_names()  # type: ignore
