import numpy as np

from pywfn import _core

core = _core.datas.mocv  # type: ignore


def get_tmat(num: int, L: np.ndarray) -> np.ndarray:
    """获取变换矩阵

    Args:
        num: 原子序数
        L: 局部坐标系

    Returns:
        np.ndarray: 变换矩阵
    """
    return core.get_tmat(num, L)  # type: ignore


def get_tmat_sph(num: int, L: np.ndarray) -> np.ndarray:
    """获取球谐变换矩阵"""
    return core.get_tmat_sph(num, L)  # type: ignore


def get_tmat_car(num: int, L: np.ndarray) -> np.ndarray:
    """获取笛卡尔变换矩阵"""
    return core.get_tmat_car(num, L)  # type: ignore


def tmat_s(L: np.ndarray) -> np.ndarray:
    """获取s轨道的变换矩阵"""
    return core.tmat_s(L)  # type: ignore


def tmat_p(L: np.ndarray) -> np.ndarray:
    """获取p轨道的变换矩阵"""
    return core.tmat_p(L)  # type: ignore


def tmat_d(L: np.ndarray) -> np.ndarray:
    """获取d轨道的变换矩阵"""
    return core.tmat_d(L)  # type: ignore


def tmat_f(L: np.ndarray) -> np.ndarray:
    """获取f轨道的变换矩阵"""
    return core.tmat_f(L)  # type: ignore


def tmat_g(L: np.ndarray) -> np.ndarray:
    """获取g轨道的变换矩阵"""
    return core.tmat_g(L)  # type: ignore


def get_nmat(lmns: list[list[int]]) -> np.ndarray:
    """获取n矩阵"""
    return core.get_nmat(lmns)  # type: ignore


def nmat_d() -> np.ndarray:
    """获取d轨道的n矩阵"""
    return core.nmat_d()  # type: ignore


def nmat_f() -> np.ndarray:
    """获取f轨道的n矩阵"""
    return core.nmat_f()  # type: ignore


def nmat_g() -> np.ndarray:
    """获取g轨道的n矩阵"""
    return core.nmat_g()  # type: ignore
