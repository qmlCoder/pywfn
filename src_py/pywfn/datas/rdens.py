"""
原子的径向电子密度

对应 rswfn 的 `datas::rdens` 模块。
"""

from pywfn import _core


def get_rad_dens_v1(atomic: int, dist: float) -> float:
    """获取指定原子的径向电子密度（方法一）

    Args:
        atomic: 原子序数
        dist: 距离原子核的距离

    Returns:
        float: 径向电子密度
    """
    return _core.datas.rdens.get_rad_dens_v1(atomic, dist)  # type: ignore


def get_rad_dens_v2(atomic: int, dist: float) -> float:
    """获取指定原子的径向电子密度（方法二）

    Args:
        atomic: 原子序数
        dist: 距离原子核的距离

    Returns:
        float: 径向电子密度
    """
    return _core.datas.rdens.get_rad_dens_v2(atomic, dist)  # type: ignore
