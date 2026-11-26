"""
挑选pi轨道
"""

import numpy as np


def judge(Ca: np.ndarray, Cb: np.ndarray, As: float) -> int:
    """判断一个轨道是否为pi轨道

    Args:
        Ca: 轨道系数A
        Cb: 轨道系数B
        As: 归一化因子

    Returns:
        int: 0表示是pi轨道
    """
    # 第一步，计算轨道的贡献
    Sa = np.sum(Ca[0] ** 2) / As
    Sb = np.sum(Cb[0] ** 2) / As
    if Sa > 0.05 or Sb > 0.05:
        return 0
    return 1
