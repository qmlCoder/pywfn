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
    def idx(self):
        return self.core.idx

    @property
    def sym(self):
        return self.core.sym

    @property
    def rad(self):
        return self.core.rad

    @property
    def row(self):
        return self.core.row

    @property
    def col(self):
        return self.core.col


ELEMS = [ELEM(core) for core in _core.datas.elems.ELEMS()]  # type: ignore


@lru_cache
def get_by_idx(idx: int) -> ELEM:
    core = _core.datas.elems.get_by_idx(idx)  # type: ignore
    return ELEM(core)


@lru_cache
def get_by_sym(sym: str) -> ELEM:
    core = _core.datas.elems.get_by_sym(sym)  # type: ignore
    return ELEM(core)
