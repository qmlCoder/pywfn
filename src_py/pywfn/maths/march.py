import numpy as np
from numpy.typing import NDArray

from pywfn import _core


def march_cube_v1(
    shape: list[int],
    grids: NDArray[np.float64],
    values: NDArray[np.float64],
    isov: float,
) -> tuple[NDArray[np.float64], NDArray[np.int32]]:
    """使用Marching Cubes算法提取等值面

    Args:
        shape: 格点形状 [nx, ny, nz]
        grids: 格点坐标
        values: 格点数值
        isov: 等值面阈值

    Returns:
        tuple: (顶点坐标, 三角形类型)
    """
    verts, types = _core.maths.march.march_cube_v1(shape, grids, values, isov)  # type: ignore
    return verts, types
