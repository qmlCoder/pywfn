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

from dataclasses import dataclass
from pywfn import utils
import numpy as np

printer = utils.Printer()

from functools import lru_cache
from collections import defaultdict
from pywfn import base

ANG2LMN = {
    0: [[0, 0, 0]],
    1: [[1, 0, 0], 
        [0, 1, 0], 
        [0, 0, 1]],
    2: [[2, 0, 0], 
        [0, 2, 0], 
        [0, 0, 2], 
        [1, 1, 0], 
        [1, 0, 1], 
        [0, 1, 1]],
}
SYM2LMN={
    "S" : [0, 0, 0],
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
LMN2SYM={
    (0, 0, 0): "S"  ,
    (1, 0, 0): "PX" ,
    (0, 1, 0): "PY" ,
    (0, 0, 1): "PZ" ,
    (2, 0, 0): "XX" ,
    (0, 2, 0): "YY" ,
    (0, 0, 2): "ZZ" ,
    (1, 1, 0): "XY" ,
    (1, 0, 1): "XZ" ,
    (0, 1, 1): "YZ" ,
}


# @dataclass
class BasisData:
    def __init__(self,atm:int,shl:int,ang:int,coe:float,alp:float) -> None:
        self.atm = atm # 原子序号，从1开始
        self.shl = shl
        self.ang = ang
        self.coe = coe
        self.alp = alp

    def __iter__(self):
        data = self.atm, self.shl, self.ang, self.alp, self.coe
        return iter(data)
    
    def __repr__(self) -> str:
        resStr=''
        # resStr+=f'{"atm":>5}{"shl":>5}{"ang":>5}{"coe":>12}{"alp":>12}\n'
        resStr+=f'{self.atm:>5}{self.shl:>5}{self.ang:>5}{self.coe:>12.4e}{self.alp:>12.4e}'
        return resStr


class Basis:
    """
    所有提前准备的基组数据
    """

    def __init__(self,mol:"base.Mol") -> None:
        """根据基组名实例化基组信息"""
        self.mol  = mol
        self.name = ''
        self.data: list[BasisData] = []
    
    @property
    def dict(self):
        basDict=defaultdict(list)
        for each in self.data:
            key=f'{each.atm}-{each.shl}-{each.ang}'
            basDict[key].append([each.alp,each.coe])
        return basDict

    # def atomics(self):
    #     atmics = [basis.sym for basis in self.data]
    #     return list(set(atmics))

    def setData(self, data: list[BasisData]):
        """设置数据
        元素,层数,角动量,指数,系数
        idx,shl,ang,exp,coe
        """
        self.data = data

    @lru_cache
    def get(self, sym: str, shl: int|None = None, ang: int|None = None) -> list[BasisData]:
        """根据原子序号获得基组"""
        basis=[]
        for each in self.data:
            s1 = self.mol.atom(each.atm).sym == sym # 第一个条件
            s2 = shl is None or shl == each.shl # 第二个条件
            s3 = ang is None or ang == each.ang # 第三个条件
            if s1 and s2 and s3:
                basis.append(each)
        assert len(basis) > 0, "没有基组信息"
        return basis

    def lmn2sym(self, lmn):
        l,m,n = lmn
        return LMN2SYM[(l,m,n)]

    def sym2ang(self, sym):
        return sum(SYM2LMN[sym])
    
    def sym2lmn(self,sym):
        return SYM2LMN[sym]

    def matMap(self): # 映射到系数矩阵的数据类型
        atos=[]
        coes=[]
        alps=[]
        lmns=[]
        for each in self.data:
            atm=each.atm
            shl=each.shl
            ang=each.ang
            coe=each.coe
            alp=each.alp
            
            for i,lmn in enumerate(ANG2LMN[ang]):
                sym=self.lmn2sym(lmn)
                key=f'{atm}-{shl}{sym}'
                ato=self.mol.atoKeys.index(key)
                atos.append(ato)
                coes.append(coe)
                alps.append(alp)
                lmns.append(lmn)
        atos=np.array(atos,dtype=np.int32)
        coes=np.array(coes,dtype=np.float64)
        alps=np.array(alps,dtype=np.float64)
        lmns=np.array(lmns,dtype=np.int32)
        idxs=np.argsort(atos)
        atos=atos[idxs].copy()
        coes=coes[idxs].copy()
        alps=alps[idxs].copy()
        lmns=lmns[idxs].copy()
        return atos,coes,alps,lmns

    def __repr__(self) -> str:
        resStr=''
        resStr+=f'{"idx":>5}{"atm":>5}{"shl":>5}{"ang":>5}{"coe":>12}{"alp":>12}\n'
        for i,each in enumerate(self.data):
            resStr+=f'{i+1:>5}{each.atm:>5}{each.shl:>5}{each.ang:>5}{each.coe:>12.4e}{each.alp:>12.4e}\n'
        return resStr
            
        