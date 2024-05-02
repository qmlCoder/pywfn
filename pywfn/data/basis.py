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
printer=utils.Printer()

from functools import lru_cache

@dataclass
class BasisData:
    atm:int # 元素
    shl:int # 壳层
    ang:int # 角动量
    exp:float # 指数
    coe:float # 系数

    def __iter__(self):
        data=self.atm,self.shl,self.ang,self.exp,self.coe
        return iter(data)

class Basis:
    """
    所有提前准备的基组数据
    """
    def __init__(self,name:str) -> None:
        """根据基组名实例化基组信息"""
        self.name=name
        self.data:list[BasisData]=None
        
    def setData(self,data:list[BasisData]):
        """设置数据
        元素,层数,角动量,指数,系数
        idx,shl,ang,exp,coe
        """
        self.data=data
    
    @lru_cache
    def get(self,atm:int,shl:int=None,ang:int=None)->list[BasisData]:
        """根据原子序号获得基组"""
        if shl is None:
            return [b for b in self.data if b.atm==atm]
        else:
            return [b for b in self.data if (b.atm==atm and b.shl==shl and b.ang==ang)]
    
    @lru_cache
    def lmn(self,ang:int)->list[list[int]]:
        """
        根据角动量获取l,m,n
        基函数中包含 x^l.y^m.z^n
        注意，角动量要与计算时的系数对应
        """
        lmns={
            0:[[0, 0, 0]],
            1:[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            2:[[2, 0, 0], [0, 2, 0], [0, 0, 2], [1, 1, 0], [1, 0, 1], [0, 1, 1]]
        }
        return lmns[ang]

    def num(self,atomic:int)->int:
        """获取基组原子对应的基函数数量"""
        data=self.get(atomic)
        return sum([len(self.lmn(ang)) for ang,_,_ in data])
    
    def lName(self,l,m,n):
        key=f'{l}{m}{n}'
        names={
            '000':'S',
            '100':'PX',
            '010':'PY',
            '001':'PZ',
            '200':'XX',
            '020':'YY',
            '002':'ZZ',
            '110':'XY',
            '101':'XZ',
            '011':'YZ',
        }
        return names[key]
    
    def numAng(self,strs):
        """将角动量符号转为数值"""
        return [[int(i) for i in s] for s in strs]