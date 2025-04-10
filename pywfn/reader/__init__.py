"""
该子模块定义各种文件的读取器
所有读取器应该都是读取分子所需要的属性？
数据尽量向量化存储
分子文件需要哪些属性，则读取器需要读取哪些属性
原子类型
原子坐标
系数矩阵
重叠矩阵
"""


from pathlib import Path


# 定义所有reader的基类
# 基类中定义各种属性的读取函数，再在子类中覆盖

"""
定义所有读取类的基类
默认读取的所有信息返回为空
不同类型的文件覆盖不同的读取函数

可以为每个分子对象分配一个Reader对象
用户可以自己分配分子属性，也可以从文件中读取
当调用分子属性且未分配时，会通过reader对象读取属性，若读取失败则报错
每一个方法都需要被子例覆盖

- 原子类型
- 原子坐标
- 轨道系数，[n,n]/[n,2n]
- 分子能量
- 原子轨道所属原子
- 原子轨道角动量类型
- 分子轨道能量
- 分子轨道占据类型 占据|非占据
- 重叠矩阵
- 波函数类型 开壳层/闭壳层
"""
from pywfn import base
from pywfn import data
from pywfn import config
import numpy as np
from pathlib import Path
from functools import cached_property
import linecache
from pywfn.utils import printer
from typing import Callable
import os
import shutil

class Reader:
    def __init__(self,path:str,cache:bool=False) -> None:
        """初始化读取器父类

        Args:
            path (str): 文件路径
            cache (bool, optional): 是否缓存. Defaults to False.
        """
        self.type:str='' # 读取器类型
        self.path:str=path
        self.cache=cache
        dataFold=Path(path).parent/f'{Path(path).suffix}#{Path(path).stem}' # 数据存储路径
        self.dataFold=f'{dataFold}' # 数据存储路径
        if config.CACHE: #如果生成缓存
            if not dataFold.exists():
                os.mkdir(f'{dataFold}')
            else:
                foldTime=dataFold.stat().st_mtime # 文件夹创建的时间
                fileTime=Path(path).stat().st_mtime #文件创建的时间
                if fileTime>foldTime: # 如果文件更新了
                    shutil.rmtree(f'{dataFold}') # 删除文件夹
                    os.mkdir(f'{dataFold}') # 重新创建
                    printer.info("清空数据缓存")
            
    
    @cached_property
    def text(self)->str:
        """获取文件的所有文本，小文件可以这么干，大文件就算了"""
        return Path(self.path).read_text()

    @property
    def fname(self)->str:
        return Path(self.path).name
    
    @cached_property
    def lineNum(self)->int:
        """获取文件行数"""
        lineNum=0
        with open(self.path,'r') as file:
            for _ in file:
                lineNum+=1
        return lineNum
    
    def getline(self,idx:int,keepEnd:bool=True)->str:
        if keepEnd:
            return linecache.getline(self.path,idx+1)
        else:
            return linecache.getline(self.path,idx+1).replace('\n','')
    
    def getlines(self,idx1:int,idx2:int)->list[str]:
        lines=[]
        for i in range(idx1,idx2):
            lines.append(self.getline(i,keepEnd=False))
        return lines

    def get_energy(self)->float:
        """获取分子能量"""
        raise ValueError("未继承的函数")

    def get_nele(self)->tuple[int,int]:
        """获取分子电子数"""
        raise ValueError("未继承的函数")
    
    def get_geome(self)->"base.Geome":
        """获取分子几何信息"""
        raise ValueError("未继承的函数")
    
    def get_basis(self)->"base.Basis":
        """获取每个原子的基组数据
        原子：壳层，角动量，系数，指数
        """
        raise ValueError("未继承的函数")
    
    def get_coefs(self)->"base.Coefs":
        """获取轨道系数"""
        raise ValueError("未继承的函数")
    
    def load_fdata(self,name:str)->np.ndarray|None:
        """获取保存的数据"""
        path=Path(self.dataFold)/name
        if not path.exists():return None
        return np.load(path)
    
    def save_fdata(self,name:str,data:np.ndarray|list):
        """保存数据"""
        path=Path(self.dataFold)/name
        np.save(path,data)


from pywfn.reader.fch import FchReader
from pywfn.reader.log import LogReader
from pywfn.reader.gjf import GjfReader
from pywfn.reader.mol import MolReader
from pywfn.reader.any import AnyReader
from pywfn.reader.mod import ModReader
from pywfn.reader.sdf import SdfReader
from pywfn.reader.xyz import XyzReader
from pywfn.reader.wfn import WfnReader

supports=[".log",".out",".fch",".gjf",".mol",".molden",".sdf",".xyz"] # 支持的文件类型

def get_reader(path:str):
    """根据输入文件的类型自动判断应该使用哪个读取器"""
    suffix=Path(path).suffix
    readers={
        '.out':LogReader,
        '.log':LogReader,
        '.gjf':GjfReader,
        '.fch':FchReader,
        '.mol':MolReader,
        '.sdf':SdfReader,
        '.xyz':XyzReader,
        '.molden':ModReader,
    }
    assert suffix in readers.keys(),f'不支持的文件类型{suffix}'
    return readers[suffix](path)