"""
日志模块

对应 rswfn 的 `loger` 模块，用于初始化日志输出。
"""

from pywfn import _core


def init_log() -> None:
    """初始化日志，输出到日志文件"""
    _core.loger.init_log()  # type: ignore


def init_stderr_log() -> None:
    """初始化日志，输出到标准错误流"""
    _core.loger.init_stderr_log()  # type: ignore
