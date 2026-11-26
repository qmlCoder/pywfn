"""
基于python的波函数分析程序
"""

from pywfn import _core


def rswfn_version() -> str:
    """返回底层 rswfn 库的版本号"""
    return _core.utils.rswfn_version()


__rswfn_version__ = rswfn_version()
