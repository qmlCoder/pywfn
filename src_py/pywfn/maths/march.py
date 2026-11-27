"""
等值面提取

对应 rswfn 的 `maths::march` 模块，提供 Marching Cubes / Marching Tetrahedra /
Surface Nets 等等值面提取算法。
"""

from typing import Callable

import numpy as np
from numpy.typing import NDArray

from pywfn import _core


class Cube:
    """立方体区域，用于等值面提取"""

    def __init__(self, p0: list[float], p1: list[float]) -> None:
        self.core = _core.maths.march.Cube(p0, p1)  # type: ignore

    def calc_step(self, count: int) -> list[float]:
        """按格点数量计算步长

        Args:
            count: 每个方向的格点数量
        """
        return self.core.calc_step(count)

    def volum(self) -> float:
        """立方体体积"""
        return self.core.volum()

    def p0(self):
        """起点坐标"""
        from pywfn.base import Vec3

        return Vec3(self.core.p0())

    def p1(self):
        """终点坐标"""
        from pywfn.base import Vec3

        return Vec3(self.core.p1())

    def grids_v1(self, count: int) -> tuple[list[int], NDArray[np.float64], float]:
        """按格点数量生成格点

        Args:
            count: 每个方向的格点数量

        Returns:
            (形状, 格点坐标, 步长)
        """
        shape, grids, step = self.core.grids_v1(count)
        return shape, np.array(grids), step

    def grids_v2(self, step: float) -> tuple[list[int], NDArray[np.float64]]:
        """按步长生成格点

        Args:
            step: 步长

        Returns:
            (形状, 格点坐标)
        """
        shape, grids = self.core.grids_v2(step)
        return shape, np.array(grids)


def march_cube_v1(
    shape: list[int],
    grids: NDArray[np.float64],
    values: NDArray[np.float64],
    isov: float,
) -> tuple[NDArray[np.float64], NDArray[np.int32]]:
    """使用 Marching Cubes 算法提取等值面

    Args:
        shape: 格点形状 [nx, ny, nz]
        grids: 格点坐标
        values: 格点数值
        isov: 等值面阈值

    Returns:
        tuple: (顶点坐标, 三角形类型)
    """
    verts, types = _core.maths.march.mc.march_cube_v1(  # type: ignore
        shape, grids, values, isov
    )
    return np.array(verts), np.array(types)


def march_cube_v2(
    func: Callable[[list[float]], float],
    cube: Cube,
    isov: float,
    step: float,
    xyzs: NDArray[np.float64] | None = None,
    limit: float | None = None,
) -> tuple[NDArray[np.float64], NDArray[np.int32]]:
    """使用 Marching Cubes 算法提取等值面（直接给定函数）

    Args:
        func: 计算函数值的回调，输入坐标列表返回浮点数
        cube: 立方体区域
        isov: 等值面阈值
        step: 步长
        xyzs: 原子坐标列表，可选
        limit: 距离限制，可选
    """
    verts, types = _core.maths.march.mc.march_cube_v2(  # type: ignore
        func, cube.core, isov, step, xyzs, limit
    )
    return np.array(verts), np.array(types)


def voxel2verts(
    voxel: NDArray[np.float64], isov: float
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """提取体素中所有格点及其数值

    Args:
        voxel: 四维体素数组
        isov: 等值面阈值
    """
    verts, acts = _core.maths.march.mc.voxel2verts(voxel, isov)  # type: ignore
    return np.array(verts), np.array(acts)


def march_tetra_v1(
    shape: list[int],
    grids: NDArray[np.float64],
    values: NDArray[np.float64],
    isov: float,
) -> tuple[NDArray[np.float64], NDArray[np.int32]]:
    """使用 Marching Tetrahedra 算法提取等值面

    Args:
        shape: 格点形状 [nx, ny, nz]
        grids: 格点坐标
        values: 格点数值
        isov: 等值面阈值
    """
    verts, types = _core.maths.march.mt.march_tetra_v1(  # type: ignore
        shape, grids, values, isov
    )
    return np.array(verts), np.array(types)


def march_tetra_v2(
    func: Callable[[list[float]], float],
    cube: Cube,
    isov: float,
    step: float,
    xyzs: NDArray[np.float64] | None = None,
    limit: float | None = None,
) -> tuple[NDArray[np.float64], NDArray[np.int32]]:
    """使用 Marching Tetrahedra 算法提取等值面（直接给定函数）"""
    verts, types = _core.maths.march.mt.march_tetra_v2(  # type: ignore
        func, cube.core, isov, step, xyzs, limit
    )
    return np.array(verts), np.array(types)


def surface_nets_v1(
    shape: list[int],
    grids: NDArray[np.float64],
    values: NDArray[np.float64],
    isov: float,
) -> tuple[NDArray[np.float64], NDArray[np.int32]]:
    """使用 Surface Nets 算法提取等值面

    Args:
        shape: 格点形状 [nx, ny, nz]
        grids: 格点坐标
        values: 格点数值
        isov: 等值面阈值
    """
    verts, types = _core.maths.march.sn.surface_nets_v1(  # type: ignore
        shape, grids, values, isov
    )
    return np.array(verts), np.array(types)


def extract_region(
    func: Callable[[list[float]], float],
    cube: Cube,
    isov: float,
    step: float,
    xyzs: NDArray[np.float64] | None = None,
    limit: float | None = None,
    extract: Callable | None = None,
) -> tuple[NDArray[np.float64], NDArray[np.int32]]:
    """从三维标量场中提取等值面区域

    Args:
        func: 计算函数值的回调
        cube: 立方体区域
        isov: 等值面阈值
        step: 步长
        xyzs: 原子坐标列表，可选
        limit: 距离限制，可选
        extract: 提取函数，接收 (形状, 格点, 数值, 阈值) 并返回 (顶点, 类型)
    """
    verts, types = _core.maths.march.extract_region(  # type: ignore
        func, cube.core, isov, step, xyzs, limit, extract
    )
    return np.array(verts), np.array(types)


def flatten(
    verts: NDArray[np.float64], faces: NDArray[np.int64]
) -> tuple[NDArray[np.float64], NDArray[np.int64]]:
    """把顶点和三角面片展平为一维数组

    Args:
        verts: 顶点坐标列表
        faces: 三角面片索引列表
    """
    vals, idxs = _core.maths.march.flatten(verts, faces)  # type: ignore
    return np.array(vals), np.array(idxs)


def grids2voxel(
    shape: list[int], grids: NDArray[np.float64], values: NDArray[np.float64]
) -> NDArray[np.float64]:
    """把格点和数值整理成四维体素数组

    Args:
        shape: 格点形状 [nx, ny, nz]
        grids: 格点坐标列表
        values: 格点数值列表
    """
    return _core.maths.march.grids2voxel(shape, grids, values)  # type: ignore


def merge_verts(
    old_verts: NDArray[np.float64],
) -> tuple[NDArray[np.float64], NDArray[np.int64]]:
    """合并重复顶点

    Args:
        old_verts: 原始顶点坐标列表
    """
    verts, faces = _core.maths.march.merge_verts(old_verts)  # type: ignore
    return np.array(verts), np.array(faces)


def smooth_isosurface(
    verts: NDArray[np.float64],
    types: NDArray[np.int32],
    isov: float,
    func: Callable[[list[float]], float],
    iterations: int,
) -> NDArray[np.float64]:
    """平滑等值面顶点

    Args:
        verts: 顶点坐标列表
        types: 顶点类型列表
        isov: 等值面阈值
        func: 计算函数值的回调
        iterations: 迭代次数
    """
    return np.array(
        _core.maths.march.smooth_isosurface(verts, types, isov, func, iterations)  # type: ignore
    )


def snap_to_isosurface(
    verts: NDArray[np.float64],
    types: NDArray[np.int32],
    isov: float,
    func: Callable[[list[float]], float],
    h: float,
) -> NDArray[np.float64]:
    """把顶点吸附到等值面上

    Args:
        verts: 顶点坐标列表
        types: 顶点类型列表
        isov: 等值面阈值
        func: 计算函数值的回调
        h: 步长
    """
    return np.array(
        _core.maths.march.snap_to_isosurface(verts, types, isov, func, h)  # type: ignore
    )


def value_merge(
    pverts: NDArray[np.float64],
    pfaces: NDArray[np.int64],
    nverts: NDArray[np.float64],
    nfaces: NDArray[np.int64],
) -> tuple[NDArray[np.float64], NDArray[np.int64], NDArray[np.int32]]:
    """按数值合并两组网格

    Args:
        pverts: 第一组顶点
        pfaces: 第一组面片
        nverts: 第二组顶点
        nfaces: 第二组面片
    """
    verts, faces, types = _core.maths.march.value_merge(  # type: ignore
        pverts, pfaces, nverts, nfaces
    )
    return np.array(verts), np.array(faces), np.array(types)


def voxel2cubes(voxel: NDArray[np.float64], isov: float) -> list[Cube]:
    """从体素数组中提取所有包含等值面的立方体

    Args:
        voxel: 四维体素数组
        isov: 等值面阈值
    """
    cubes: list[Cube] = []
    for core in _core.maths.march.voxel2cubes(voxel, isov):  # type: ignore
        cube = Cube.__new__(Cube)
        cube.core = core
        cubes.append(cube)
    return cubes
