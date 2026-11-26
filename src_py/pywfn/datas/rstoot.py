from pywfn import _core


def rys_root1(x: float) -> tuple[list[float], list[float]]:
    """计算Rys多项式的一阶根"""
    return _core.datas.rsroot.rys_root1(x)  # type: ignore


def rys_root2(x: float) -> tuple[list[float], list[float]]:
    """计算Rys多项式的二阶根"""
    return _core.datas.rsroot.rys_root2(x)  # type: ignore


def rys_root3(x: float) -> tuple[list[float], list[float]]:
    """计算Rys多项式的三阶根"""
    return _core.datas.rsroot.rys_root3(x)  # type: ignore


def rys_root4(x: float) -> tuple[list[float], list[float]]:
    """计算Rys多项式的四阶根"""
    return _core.datas.rsroot.rys_root4(x)  # type: ignore


def rys_root5(x: float) -> tuple[list[float], list[float]]:
    """计算Rys多项式的五阶根"""
    return _core.datas.rsroot.rys_root5(x)  # type: ignore
