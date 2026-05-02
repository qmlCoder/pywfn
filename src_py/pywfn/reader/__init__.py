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
from pywfn.base.atoms import Atoms
from pywfn.base.basis import Basis
from pywfn.base.coefs import Coefs


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
from pathlib import Path

from pywfn.reader.fch import FchReader
from pywfn.reader.log import LogReader
from pywfn.reader.gjf import GjfReader
from pywfn.reader.mol import MolReader
from pywfn.reader.non import NonReader
from pywfn.reader.mde import MdeReader
from pywfn.reader.sdf import SdfReader
from pywfn.reader.xyz import XyzReader
from pywfn.reader.cub import CubReader
# from pywfn.reader.wfn import WfnReader

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
        '.molden':MdeReader,
    }
    assert suffix in readers.keys(),f'不支持的文件类型{suffix}'
    return readers[suffix](path)