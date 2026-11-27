"""
HMO（Hückel 分子轨道）积分参数

对应 rswfn 的 `datas::dhmo2` 模块。
"""

from pywfn import _core


def get_paras(Zi: int, Zj: int, Pi: int, Pj: int) -> list[list[float]]:
    """获取 HMO 积分参数

    Args:
        Zi: 原子 i 的电荷数
        Zj: 原子 j 的电荷数
        Pi: 原子 i 的主量子数
        Pj: 原子 j 的主量子数

    Returns:
        list[list[float]]: 59 组参数，每组 4 个值
    """
    return _core.datas.dhmo2.get_paras(Zi, Zj, Pi, Pj)  # type: ignore
