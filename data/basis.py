"""
一个原子的基组应该有一下层级

不同价层
    指数(一维数组)
    系数(二维数组)
不同角动量对应不同的mnl
所有可能的m+n+l=角动量决定了基函数的数量

"""

from dataclasses import dataclass
from pywfn import utils
import numpy as np

printer = utils.Printer()

from functools import lru_cache


@dataclass
class BasisData:
    atmic: int  # 元素
    shl: int  # 壳层
    ang: int  # 角动量
    exp: float  # 指数
    coe: float  # 系数
    # 对应的矩阵

    def __iter__(self):
        data = self.atmic, self.shl, self.ang, self.exp, self.coe
        return iter(data)


class Basis:
    """
    所有提前准备的基组数据
    """

    def __init__(self, name: str) -> None:
        """根据基组名实例化基组信息"""
        self.name = name
        self.data: list[BasisData] = None

    def atomics(self):
        atmics = [basis.atmic for basis in self.data]
        return list(set(atmics))

    def setData(self, data: list[BasisData]):
        """设置数据
        元素,层数,角动量,指数,系数
        idx,shl,ang,exp,coe
        """
        self.data = data

    @lru_cache
    def get(self, atmic: int, shell: int = None, ang: int = None) -> list[BasisData]:
        """根据原子序号获得基组"""
        if shell is None:
            basis = [b for b in self.data if b.atmic == atmic]
        else:
            basis = [
                b
                for b in self.data
                if (b.atmic == atmic and b.shl == shell and b.ang == ang)
            ]
        assert len(basis) > 0, "没有基组信息"
        return basis

    @staticmethod
    def ang2lmn(ang: int) -> list[list[int]]:
        """
        根据角动量获取l,m,n
        基函数中包含 x^l.y^m.z^n
        注意，角动量要与计算时的系数对应
        """
        lmns = {
            0: [[0, 0, 0]],
            1: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            2: [[2, 0, 0], [0, 2, 0], [0, 0, 2], [1, 1, 0], [1, 0, 1], [0, 1, 1]],
        }
        return lmns[ang]

    @staticmethod
    def sym2lmn(sym: str)->list[int]:
        lmnMap = {
            "S": [0, 0, 0],
            "PX": [1, 0, 0],
            "PY": [0, 1, 0],
            "PZ": [0, 0, 1],
            "XX": [2, 0, 0],
            "YY": [0, 2, 0],
            "ZZ": [0, 0, 2],
            "XY": [1, 1, 0],
            "XZ": [1, 0, 1],
            "YZ": [0, 1, 1],
        }
        return lmnMap[sym]
    
    @staticmethod
    def sym2ang(sym: str)->int:
        lmn=Basis.sym2lmn(sym)
        return sum(lmn)

    def num(self, atomic: int) -> int:
        """获取基组原子对应的基函数数量"""
        data = self.get(atomic)
        return sum([len(self.lmn(ang)) for ang, _, _ in data])

    def lmn2sym(self, lmn):
        l,m,n = lmn
        key = f"{l}{m}{n}"
        names = {
            "000": "S",
            "100": "PX",
            "010": "PY",
            "001": "PZ",
            "200": "XX",
            "020": "YY",
            "002": "ZZ",
            "110": "XY",
            "101": "XZ",
            "011": "YZ",
        }
        return names[key]

    def matMap(self,atmic:int):
        from collections import defaultdict
        alpsDict=defaultdict(list)
        coesDict=defaultdict(list)
        keys=[]
        for each in self.data:
            if each.atmic!=atmic:continue
            shl=each.shl
            lmns = self.ang2lmn(each.ang)
            for lmn in lmns:
                l,m,n=lmn
                key=f'{each.atmic},{shl},{l},{m},{n}'
                if key not in keys:keys.append(key)
                alpsDict[key].append(each.exp)
                coesDict[key].append(each.coe)

        ncgs=[len(alpsDict[key]) for key in keys] # 每一个原子轨道对应的收缩数量
        alps=[]
        coes=[]
        for i,key in enumerate(keys):
            alps.append(alpsDict[key])
            coes.append(coesDict[key])
        assert len(ncgs)==len(alps),"长度需一致"
        return ncgs,alps,coes
        