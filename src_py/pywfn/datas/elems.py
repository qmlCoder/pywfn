"""
记录元素数据
- 元素序号
- 元素符号
- 原子半径 https://chemistry-europe.onlinelibrary.wiley.com/doi/full/10.1002/chem.200800987
"""

from functools import lru_cache

from pywfn import _core


class ELEM:
    def __init__(self, core: _core.datas.elems.ELEM) -> None:  # type: ignore
        self.core = core

    @property
    def idx(self) -> int:
        """原子序数"""
        return self.core.idx

    @property
    def sym(self) -> str:
        """元素符号"""
        return self.core.sym

    @property
    def rad(self) -> float:
        """原子半径"""
        return self.core.rad

    @property
    def row(self) -> int:
        """周期数"""
        return self.core.row

    @property
    def col(self) -> int:
        """族数"""
        return self.core.col

    @property
    def color(self) -> str:
        """颜色"""
        return self.core.color

    @property
    def oxi_num(self) -> list[int]:
        """氧化数

        `OXIS` 从 H 开始按 0 下标排列，而 `idx` 是 1 基的原子序数，因此需要减一；
        对 `X`(-1)、`Bq`(0) 这类伪原子没有氧化数，返回空列表。
        """
        if self.core.idx < 1:
            return []
        return OXIS[self.core.idx - 1]


ELEMS = [ELEM(core) for core in _core.datas.elems.ELEMS()]  # type: ignore


@lru_cache
def get_by_idx(idx: int) -> ELEM:
    core = _core.datas.elems.get_by_idx(idx)  # type: ignore
    return ELEM(core)


@lru_cache
def get_by_sym(sym: str) -> ELEM:
    core = _core.datas.elems.get_by_sym(sym)  # type: ignore
    return ELEM(core)


@lru_cache
def NLMS() -> list[tuple[int, int, int]]:
    core = _core.datas.elems.NLMS()  # type: ignore
    return core


@lru_cache
def ELES() -> list[list[int]]:
    core = _core.datas.elems.ELES()  # type: ignore
    return core


def get_ab_nele(atm: int) -> list[int]:
    """获取指定元素的α电子数和β电子数

    参数:
        atm: 原子序号

    返回:
        [nela, nelb] 分别为α电子数和β电子数
    """
    return _core.datas.elems.get_ab_nele(atm)  # type: ignore


OXIS: list[list[int]] = _core.datas.elems.OXIS()  # type: ignore


@lru_cache
def sym2atm(sym: str) -> int:
    """根据元素符号返回原子序数"""
    return _core.datas.elems.sym2atm(sym)  # type: ignore


@lru_cache
def atm2sym(atm: int) -> str:
    """根据原子序数返回元素符号"""
    return _core.datas.elems.atm2sym(atm)  # type: ignore


@lru_cache
def zeff(atm: int, neb: int) -> float:
    """返回指定原子的有效核电荷数（用于pi电子分析）

    Args:
        atm: 原子序号
        neb: 参与共轭的电子数
    """
    return _core.datas.elems.zeff(atm, neb)  # type: ignore


@lru_cache
def atm2rad(idx: int) -> float:
    """根据原子序数返回原子半径"""
    return _core.datas.elems.atm2rad(idx)  # type: ignore


@lru_cache
def atm2row(atm: int) -> int:
    """根据原子序数返回周期数"""
    return _core.datas.elems.atm2row(atm)  # type: ignore


def ELEM_DATA() -> list[tuple[int, str, float, int, int, str]]:
    """返回原始的元组形式元素数据表

    Returns:
        (原子序数, 元素符号, 原子半径, 周期数, 族数, 颜色) 的列表
    """
    return _core.datas.elems.ELEM_DATA()  # type: ignore

def get_by_sym_no_cache(sym: str) -> ELEM:
    """通过元素符号查找元素信息（不使用缓存）"""
    return ELEM(_core.datas.elems.get_by_sym_no_cache(sym))  # type: ignore


def get_by_sym_prime_cache(sym: str) -> ELEM:
    """通过元素符号查找元素信息（预填充缓存）"""
    return ELEM(_core.datas.elems.get_by_sym_prime_cache(sym))  # type: ignore


def get_by_idx_no_cache(idx: int) -> ELEM:
    """按原子序数查找元素信息（不使用缓存）"""
    return ELEM(_core.datas.elems.get_by_idx_no_cache(idx))  # type: ignore


def get_by_idx_prime_cache(idx: int) -> ELEM:
    """按原子序数查找元素信息（预填充缓存）"""
    return ELEM(_core.datas.elems.get_by_idx_prime_cache(idx))  # type: ignore


def atm2sym_no_cache(atm: int) -> str:
    """根据原子序数返回元素符号（不使用缓存）"""
    return _core.datas.elems.atm2sym_no_cache(atm)  # type: ignore


def atm2sym_prime_cache(atm: int) -> str:
    """根据原子序数返回元素符号（预填充缓存）"""
    return _core.datas.elems.atm2sym_prime_cache(atm)  # type: ignore


def atm2rad_no_cache(idx: int) -> float:
    """根据原子序数返回原子半径（不使用缓存）"""
    return _core.datas.elems.atm2rad_no_cache(idx)  # type: ignore


def atm2rad_prime_cache(idx: int) -> float:
    """根据原子序数返回原子半径（预填充缓存）"""
    return _core.datas.elems.atm2rad_prime_cache(idx)  # type: ignore


def atm2row_no_cache(atm: int) -> int:
    """根据原子序数返回周期数（不使用缓存）"""
    return _core.datas.elems.atm2row_no_cache(atm)  # type: ignore


def atm2row_prime_cache(atm: int) -> int:
    """根据原子序数返回周期数（预填充缓存）"""
    return _core.datas.elems.atm2row_prime_cache(atm)  # type: ignore


def zeff_no_cache(atm: int, neb: int) -> float:
    """返回指定原子的有效核电荷数（不使用缓存）"""
    return _core.datas.elems.zeff_no_cache(atm, neb)  # type: ignore


def zeff_prime_cache(atm: int, neb: int) -> float:
    """返回指定原子的有效核电荷数（预填充缓存）"""
    return _core.datas.elems.zeff_prime_cache(atm, neb)  # type: ignore