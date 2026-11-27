"""
函数优化

对应 rswfn 的 `maths::optim` 模块，提供共轭梯度、牛顿法等多种优化算法。
"""

from typing import Callable, Union

import numpy as np
from numpy.typing import NDArray

from pywfn import _core


class FuncEval:
    """一次函数求值的结果"""

    def __init__(
        self,
        value: float,
        grad: NDArray[np.float64] | None = None,
        hessian: NDArray[np.float64] | None = None,
    ) -> None:
        self.core = _core.maths.optim.FuncEval(value, grad, hessian)  # type: ignore

    @staticmethod
    def from_value(value: float) -> "FuncEval":
        """仅给定函数值"""
        feval = FuncEval.__new__(FuncEval)
        feval.core = _core.maths.optim.FuncEval.from_value(value)  # type: ignore
        return feval

    @staticmethod
    def from_grad(value: float, grad: NDArray[np.float64]) -> "FuncEval":
        """给定函数值和梯度"""
        feval = FuncEval.__new__(FuncEval)
        feval.core = _core.maths.optim.FuncEval.from_grad(value, grad)  # type: ignore
        return feval

    @staticmethod
    def from_hessian(value: float, hessian: NDArray[np.float64]) -> "FuncEval":
        """给定函数值和 Hessian 矩阵"""
        feval = FuncEval.__new__(FuncEval)
        feval.core = _core.maths.optim.FuncEval.from_hessian(value, hessian)  # type: ignore
        return feval

    @property
    def value(self) -> float:
        """函数值"""
        return self.core.value

    @property
    def grad(self) -> list[float] | None:
        """梯度"""
        return self.core.grad

    @property
    def hessian(self) -> NDArray[np.float64] | None:
        """Hessian 矩阵"""
        return self.core.hessian


OptCallable = Callable[[list[float]], Union[float, FuncEval, tuple]]


def minimize(
    func: OptCallable,
    paras: list[float],
    conv_thresh: float = 1e-6,
    method: str = "cg",
) -> list[float]:
    """最小化目标函数

    Args:
        func: 回调函数，接收参数列表，返回函数值（或 FuncEval、(value, grad, hessian) 元组）
        paras: 初始参数列表
        conv_thresh: 收敛阈值
        method: 优化算法，可选 "cg"、"dn"、"nt"、"qn"、"sd"

    Returns:
        list[float]: 优化后的参数列表
    """
    calers = {
        "cg": _core.maths.optim.cg.minimize,  # type: ignore
        "dn": _core.maths.optim.dn.minimize,  # type: ignore
        "nt": _core.maths.optim.nt.minimize,  # type: ignore
        "qn": _core.maths.optim.qn.minimize,  # type: ignore
        "sd": _core.maths.optim.sd.minimize,  # type: ignore
    }
    if method not in calers:
        raise ValueError(f"未知的优化算法: {method}")
    return calers[method](func, paras, conv_thresh)


def cg(func: OptCallable, paras: list[float], conv_thresh: float = 1e-6) -> list[float]:
    """使用共轭梯度法最小化目标函数"""
    return _core.maths.optim.cg.minimize(func, paras, conv_thresh)  # type: ignore


def dn(func: OptCallable, paras: list[float], conv_thresh: float = 1e-6) -> list[float]:
    """使用阻尼牛顿法最小化目标函数"""
    return _core.maths.optim.dn.minimize(func, paras, conv_thresh)  # type: ignore


def nt(func: OptCallable, paras: list[float], conv_thresh: float = 1e-6) -> list[float]:
    """使用牛顿法最小化目标函数"""
    return _core.maths.optim.nt.minimize(func, paras, conv_thresh)  # type: ignore


def qn(func: OptCallable, paras: list[float], conv_thresh: float = 1e-6) -> list[float]:
    """使用拟牛顿法最小化目标函数"""
    return _core.maths.optim.qn.minimize(func, paras, conv_thresh)  # type: ignore


def sd(func: OptCallable, paras: list[float], conv_thresh: float = 1e-6) -> list[float]:
    """使用最速下降法最小化目标函数"""
    return _core.maths.optim.sd.minimize(func, paras, conv_thresh)  # type: ignore
