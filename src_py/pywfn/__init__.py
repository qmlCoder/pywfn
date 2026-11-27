"""
基于python的波函数分析程序
"""

from pywfn import _core


def rswfn_version() -> str:
    """返回底层 rswfn 库的版本号"""
    return _core.utils.rswfn_version()


__rswfn_version__ = rswfn_version()

# 按 rswfn 的模块层级导出各个子模块，方便直接通过 `pywfn.<模块>` 访问
from pywfn import (  # noqa: E402
    base,
    datas,
    editor,
    loger,
    maths,
    props,
    qmcal,
    reader,
    rsbse,
    rsint,
    rsxcs,
    space,
    tools,
    utils,
    writer,
)

__all__ = [
    "base",
    "datas",
    "editor",
    "loger",
    "maths",
    "props",
    "qmcal",
    "reader",
    "rsbse",
    "rsint",
    "rsxcs",
    "space",
    "tools",
    "utils",
    "writer",
]
