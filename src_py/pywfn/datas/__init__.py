"""
内置数据

数据表统一来自底层 rswfn 的 `datas` 模块，避免在 Python 层重复维护一份。
"""

from pywfn import _core
from pywfn.datas.integ import RHM as Rhm
from pywfn.datas.integ import WHM as Whm

start = """
欢迎使用 pywfn(1.0.7): 基于python的波函数分析工具
文档 https://www.xiaofei911.top/mkdocs/pywfn/
按 q 即可退出程序，按 Enter 可返回上一级
按照提示输入编号即可轻松使用本程序 (●ˇ∀ˇ●)
"""

# 阶乘相关的常数表
FACS: list[float] = _core.datas.FACS()  # type: ignore
