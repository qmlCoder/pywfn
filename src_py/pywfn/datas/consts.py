"""
定义一些物理常数

对应 rswfn 的 `datas::consts` 模块。
"""

from pywfn import _core

R = 8.314  # 理想气体常数

# 玻尔半径（单位：埃），取自底层 rswfn 的常数
BOHR = _core.datas.consts.BOHR()  # type: ignore

# 单位转换
