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

from pywfn import _core


# class MapAto:
#     def __init__(self, core:_core.base.MapAto) -> None:  # type: ignore
#         self.core = core  # type: ignore

#     def len(self):
#         return self.core.len()

#     def coes_norm(self):
#         return self.core.coes_norm()


# class BasisData:
#     def __init__(self, core:_core.base.BasisData) -> None: # type: ignore
#         self.core = core  # type: ignore


class Basis:
    """
    存储基组数据的类
    """

    def __init__(self, core: _core.base.Basis) -> None:  # type: ignore
        self.core = core

    # @staticmethod
    # def new(data: list[BasisData]) -> "Basis":
    #     _core = _core.base.Basis(data)  # type: ignore
    #     return Basis(_core)

    def build(self, data):
        self.core.build(data)

    def __repr__(self) -> str:
        return self.core.__repr__()
