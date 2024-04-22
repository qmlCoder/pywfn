from functools import lru_cache,cached_property

from itertools import product
from pywfn.data.basis import Basis
from pywfn.data.elements import Elements


start=\
"""
欢迎使用 pywfn : 基于python的波函数分析工具
文档 https://www.xiaofei911.top/mkdocs/pywfn/
按q 即可退出程序，如无特殊提示，按 Enter 可返回上一级
按照提示输入编号即可轻松使用本程序
"""