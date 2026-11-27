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


class Ato:
    """单个原子轨道（基函数）"""

    def __init__(self, core: _core.base.Ato) -> None:  # type: ignore
        self.core = core

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
    def sym(self) -> str:
        """轨道符号"""
        return self.core.sym

    @property
    def loc(self) -> int:
        """起始位置"""
        return self.core.loc

    @property
    def form(self) -> str:
        """基函数形式"""
        return self.core.form

    @property
    def alps(self) -> list[float]:
        """高斯指数列表"""
        return self.core.alps

    @property
    def coes(self) -> list[float]:
        """收缩系数列表"""
        return self.core.coes

    def coes_norm(self) -> list[float]:
        """归一化后的收缩系数"""
        return self.core.coes_norm()

    def lmn(self) -> list[int]:
        """量子数 (l, m, n)"""
        return self.core.lmn()


class Bas:
    """单个高斯基函数（primitive）的原始数据"""

    def __init__(self, core: _core.base.Bas) -> None:  # type: ignore
        self.core = core

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

    def locs(self, form: str = "raw") -> list[int]:
        """每个基函数的起始位置"""
        return self.core.locs(form)


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

    def coes_norm(self) -> list[float]:
        """归一化后的收缩系数列表"""
        return self.core.coes_norm()


class Shls:
    def __init__(self, core: _core.base.Shls) -> None:  # type: ignore
        self.core = core

    def get(self, i: int) -> Shl:
        return Shl(self.core.get(i))

    def len(self) -> int:
        return self.core.len()

    def uplo(self, form: str = "raw") -> list[list[int]]:
        return self.core.uplo(form)

    def nato(self, form: str = "raw") -> int:
        """指定基函数形式的基函数数量"""
        return self.core.nato(form)

    def atm_uplo(self, form: str = "raw") -> list[list[int]]:
        """指定基函数形式的原子对索引列表"""
        return self.core.atm_uplo(form)


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

    def build(self, shls: list[Shl]) -> None:
        """重建基组

        Args:
            shls: 壳层列表
        """
        self.core.build([shl.core for shl in shls])

    def is_none(self) -> bool:
        """基组是否为空"""
        return self.core.is_none()

    def int_form(self) -> str:
        """积分使用的基函数形式"""
        return self.core.int_form()

    def unify_form(self, form: str) -> "Basis":
        """把基组统一到指定基函数形式

        Args:
            form: 目标基函数形式

        Returns:
            Basis: 新的基组对象
        """
        return Basis(self.core.unify_form(form), self.mole)

    def get_by_atm(self, atm: int) -> list[Shl]:
        """获取指定原子的所有壳层

        Args:
            atm: 原子索引
        """
        return [Shl(core) for core in self.core.get_by_atm(atm)]

    def ato_keys(self) -> list[str]:
        """获取所有原子轨道的键"""
        return self.core.ato_keys()

    def nbas(self) -> int:
        """基函数数量"""
        return self.core.nbas()

    def nshl(self) -> int:
        """壳层数量"""
        return self.core.nshl()

    def equal(self, other: "Basis") -> bool:
        """判断两个基组是否相同

        Args:
            other: 另一个基组
        """
        return self.core.equal(other.core)

    def __repr__(self) -> str:
        return self.core.__repr__()
