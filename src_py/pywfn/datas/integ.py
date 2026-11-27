"""
积分相关的数据表

对应 rswfn 的 `datas::integ` 模块，提供径向哈密顿矩阵和权重矩阵。
"""

import numpy as np

from pywfn import _core

# 径向哈密顿矩阵
RHM: np.ndarray = _core.datas.integ.RHM()  # type: ignore

# 径向权重矩阵
WHM: np.ndarray = _core.datas.integ.WHM()  # type: ignore
