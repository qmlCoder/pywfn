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
from pywfn.reader.lutils import Reader
from pywfn.reader.fch import FchReader
from pywfn.reader.log import LogReader
from pywfn.reader.gjf import GjfReader
from pywfn.reader.mol import MolReader

from pathlib import Path


# 定义所有reader的基类
# 基类中定义各种属性的读取函数，再在子类中覆盖



def get_reader(path:str):
    """根据输入文件的类型自动判断应该使用哪个读取器"""
    suffix=Path(path).suffix
    readers={
        '.out':LogReader,
        '.log':LogReader,
        '.gjf':GjfReader
    }
    assert suffix in readers.keys(),'不支持的读取'
    return readers[suffix](path)
    