"""
定义所有读取类的基类
默认读取的所有信息返回为空
不同类型的文件覆盖不同的读取函数

可以为每个分子对象分配一个Reader对象
用户可以自己分配分子属性，也可以从文件中读取
当调用分子属性且未分配时，会通过reader对象读取属性，若读取失败则报错
每一个方法都需要被子例覆盖
"""
from pywfn import base
from pywfn import data
import numpy as np
from pathlib import Path

class Reader:
    def __init__(self,path) -> None:
        self.path:str=path
        self.text=Path(self.path).read_text(encoding='utf-8')
        self.lines=self.text.splitlines(keepends=False)

    def get_coords(self)->np.ndarray:
        """原子坐标[n,3]"""
        raise

    def get_symbols(self)->list[str]:
        """原子符号[n]"""
        raise

    def get_CM(self)->np.ndarray:
        """系数矩阵"""
        raise

    def get_energy(self)->float:
        raise

    def get_obtEngs(self)->list[float]:
        """获取分子轨道能量"""
        raise

    def get_obtTypes(self)->list[str]:
        """获取轨道类型"""
        raise

    def get_obtLayer(self)->list[str]:
        """获取原子轨道层类型"""
        raise

    def get_obtAtoms(self)->list[int]:
        """获取轨道系数每一行对应的原子"""
        raise

    def get_SM(self)->np.ndarray:
        """获取重叠矩阵"""
        raise

    def get_charge(self)->int:
        """获取分子电荷"""
        raise

    def get_spin(self)->int:
        """获取分子自旋"""
        raise
    
    def get_basis(self)->"data.Basis":
        """获取基组数据[n,4]
        元素,层数,角动量,指数,系数
        """
        raise