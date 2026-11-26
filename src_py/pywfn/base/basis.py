"""
一个原子的基组应该有一下层级
真实存储数据的时候只需要针对每种元素存储就可以了

不同价层
    指数(一维数组)
    系数(二维数组)
不同角动量对应不同的mnl
所有可能的m+n+l=角动量决定了基函数的数量

分子中相同元素的两个原子使用的基函数可能不同，因此要存储每个原子的数据
"""

from functools import cached_property

import pywfn
from pywfn import _core


class Atos:
    def __init__(self, _core: _core.base.Atos) -> None:  # type: ignore
        self.core = _core

    def len(self, form: str = "raw") -> int:
        """基函数数量"""
        return self.core.len(form)

    def atms(self, form: str = "raw") -> list[int]:
        """每个基函数对应的原子索引"""
        return self.core.atms(form)

    def shls(self, form: str = "raw") -> list[int]:
        """每个基函数对应的壳层索引"""
        return self.core.shls(form)

    def syms(self, form: str = "raw") -> list[str]:
        """每个基函数的符号"""
        return self.core.syms(form)


class Shl:
    def __init__(self, core: _core.base.Shl) -> None:  # type: ignore
        self.core = core

    @cached_property
    def atm(self) -> int:
        return self.core.atm

    @cached_property
    def shl(self) -> int:
        return self.core.shl

    @cached_property
    def ang(self) -> int:
        return self.core.ang

    @cached_property
    def loc(self) -> int:
        return self.core.loc

    @cached_property
    def form(self) -> str:
        return self.core.form

    @cached_property
    def alps(self) -> list[float]:
        return self.core.alps

    @cached_property
    def coes(self) -> list[float]:
        return self.core.coes


class Shls:
    def __init__(self, core: _core.base.Shls) -> None:  # type: ignore
        self.core = core

    def get(self, i: int) -> Shl:
        return Shl(self.core.get(i))

    def len(self) -> int:
        return self.core.len()

    def uplo(self, form: str = "raw") -> list[list[int]]:
        return self.core.uplo(form)


class Basis:
    """
    存储基组数据的类
    """

    def __init__(self, core: _core.base.Basis, mole: "pywfn.base.Mole") -> None:  # type: ignore
        self.core = core
        self.mole = mole

    @property
    def atos(self):
        return Atos(self.core.atos())

    def get_smat(self, form: str = "raw"):
        return self.core.get_smat(form, self.mole.core)

    def atm_uplo(self, form: str = "raw") -> list[list[int]]:
        return self.core.atm_uplo(form)

    @cached_property
    def shls(self) -> Shls:
        return Shls(self.core.shls())

    def __repr__(self) -> str:
        return self.core.__repr__()
