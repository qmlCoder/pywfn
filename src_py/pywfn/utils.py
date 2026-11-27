"""
通用工具函数

对应 rswfn 的 `utils` 模块，提供文件名/文件类型解析、矩阵比较与打印等常用功能。
"""

import numpy as np
from numpy.typing import NDArray

from pywfn import _core


def rswfn_version() -> str:
    """返回底层 rswfn 库的版本号

    Returns:
        str: rswfn 的版本号
    """
    return _core.utils.rswfn_version()


def all_same(vals0: list[float], vals1: list[float]) -> bool:
    """逐个比较两个列表对应位置的值是否完全相同

    Args:
        vals0: 第一个列表
        vals1: 第二个列表

    Returns:
        bool: 长度相同且每个元素都相等时返回 True
    """
    return _core.utils.all_same(vals0, vals1)  # type: ignore


def comp_matrix(mat1: NDArray[np.float64], mat2: NDArray[np.float64], show: int = 0) -> None:
    """比较两个矩阵，打印首个不相同的元素位置

    Args:
        mat1: 第一个矩阵
        mat2: 第二个矩阵
        show: 最大打印个数，非正数表示不打印
    """
    _core.utils.comp_matrix(mat1, mat2, show)  # type: ignore


def file_name(path: str) -> str:
    """获取文件路径中的文件名（不含目录部分）

    Args:
        path: 文件路径

    Returns:
        str: 文件名
    """
    return _core.utils.file_name(path)  # type: ignore


def file_type(path: str) -> str:
    """获取文件路径的扩展名

    Args:
        path: 文件路径

    Returns:
        str: 文件扩展名
    """
    return _core.utils.file_type(path)  # type: ignore


def fort_form(val: float, pec: int) -> str:
    """按 Fortran 风格格式化浮点数

    Args:
        val: 需要格式化的值
        pec: 小数位数

    Returns:
        str: 格式化后的字符串
    """
    return _core.utils.fort_form(val, pec)  # type: ignore


def has_nan(matrix: NDArray[np.float64]) -> bool:
    """判断矩阵中是否存在 NaN

    Args:
        matrix: 待检查的矩阵

    Returns:
        bool: 是否存在 NaN
    """
    return _core.utils.has_nan(matrix)  # type: ignore


def posi_matrix(mat: NDArray[np.float64]) -> bool:
    """判断矩阵是否所有元素都为正

    Args:
        mat: 待检查的矩阵

    Returns:
        bool: 是否所有元素都为正
    """
    return _core.utils.posi_matrix(mat)  # type: ignore


def same_matrix(mat1: NDArray[np.float64], mat2: NDArray[np.float64]) -> bool:
    """判断两个矩阵是否完全相同

    Args:
        mat1: 第一个矩阵
        mat2: 第二个矩阵

    Returns:
        bool: 两个矩阵是否相同
    """
    return _core.utils.same_matrix(mat1, mat2)  # type: ignore


def show_matrix(mat: NDArray[np.float64]) -> None:
    """打印矩阵内容

    Args:
        mat: 需要打印的矩阵
    """
    _core.utils.show_matrix(mat)  # type: ignore


def symm_matrix(mat: NDArray[np.float64]) -> bool:
    """判断矩阵是否对称

    Args:
        mat: 待检查的矩阵

    Returns:
        bool: 矩阵是否对称
    """
    return _core.utils.symm_matrix(mat)  # type: ignore
