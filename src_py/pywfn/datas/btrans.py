"""
记录基函数变换矩阵
http://sobereva.com/97 感谢卢天大大
"""

import numpy as np

from pywfn import _core

DMAT: np.ndarray = _core.datas.btrans.DMAT()  # type: ignore [ 6, 5]
FMAT: np.ndarray = _core.datas.btrans.FMAT()  # type: ignore [10, 7]
GMAT: np.ndarray = _core.datas.btrans.GMAT()  # type: ignore [15, 9]
HMAT: np.ndarray = _core.datas.btrans.HMAT()  # type: ignore [21,11]

def sph2car_tmat(ang: int) -> np.ndarray:
    return _core.datas.btrans.sph2car_tmat(ang) # type: ignore

def car2sph_tmat(ang: int) -> np.ndarray:
    return _core.datas.btrans.car2sph_tmat(ang) # type: ignore

def car2sph_fmat(
    smat_car: np.ndarray,
    angs: list[int],
) -> np.ndarray:
    return _core.datas.btrans.car2sph_fmat(smat_car, angs) # type: ignore

# SYMS
MIX_S_SYMS: list[str] = _core.datas.btrans.MIX_S_SYMS()  # type: ignore
MIX_P_SYMS: list[str] = _core.datas.btrans.MIX_P_SYMS()  # type: ignore
CAR_D_SYMS: list[str] = _core.datas.btrans.CAR_D_SYMS()  # type: ignore
CAR_F_SYMS: list[str] = _core.datas.btrans.CAR_F_SYMS()  # type: ignore
CAR_G_SYMS: list[str] = _core.datas.btrans.CAR_G_SYMS()  # type: ignore
CAR_H_SYMS: list[str] = _core.datas.btrans.CAR_H_SYMS()  # type: ignore
SPH_D_SYMS: list[str] = _core.datas.btrans.SPH_D_SYMS()  # type: ignore
SPH_F_SYMS: list[str] = _core.datas.btrans.SPH_F_SYMS()  # type: ignore
SPH_G_SYMS: list[str] = _core.datas.btrans.SPH_G_SYMS()  # type: ignore
SPH_H_SYMS: list[str] = _core.datas.btrans.SPH_H_SYMS()  # type: ignore

# LMNS
MIX_S_LMNS: list[list[int]] = _core.datas.btrans.MIX_S_LMNS()  # type: ignore
MIX_P_LMNS: list[list[int]] = _core.datas.btrans.MIX_P_LMNS()  # type: ignore
CAR_D_LMNS: list[list[int]] = _core.datas.btrans.CAR_D_LMNS()  # type: ignore
CAR_F_LMNS: list[list[int]] = _core.datas.btrans.CAR_F_LMNS()  # type: ignore
CAR_G_LMNS: list[list[int]] = _core.datas.btrans.CAR_G_LMNS()  # type: ignore
CAR_H_LMNS: list[list[int]] = _core.datas.btrans.CAR_H_LMNS()  # type: ignore

SYMS: list[str] = _core.datas.btrans.SYMS()  # type: ignore
LMNS: list[list[int]] = _core.datas.btrans.LMNS()  # type: ignore


def ang2lmns(ang: int) -> list[list[int]]:
    """获取指定角动量的量子数 (l, m, n) 列表

    Args:
        ang: 角动量量子数
    """
    return _core.datas.btrans.ang2lmns(ang)  # type: ignore


def ang2syms(ang: int, form: str = "sph") -> list[str]:
    """获取指定角动量的对称性标签列表

    Args:
        ang: 角动量量子数
        form: 基函数形式（"sph" 或 "car"）
    """
    return _core.datas.btrans.ang2syms(ang, form)  # type: ignore


def get_car_syms() -> list[str]:
    """获取笛卡尔型基函数的对称性标签"""
    return _core.datas.btrans.get_car_syms()  # type: ignore


def get_car_syms_no_cache() -> list[str]:
    """获取笛卡尔型基函数的对称性标签（不使用缓存）"""
    return _core.datas.btrans.get_car_syms_no_cache()  # type: ignore


def get_car_syms_prime_cache() -> list[str]:
    """获取笛卡尔型基函数的对称性标签（预填充缓存）"""
    return _core.datas.btrans.get_car_syms_prime_cache()  # type: ignore


def get_sph_syms() -> list[str]:
    """获取球谐型基函数的对称性标签"""
    return _core.datas.btrans.get_sph_syms()  # type: ignore


def get_sph_syms_no_cache() -> list[str]:
    """获取球谐型基函数的对称性标签（不使用缓存）"""
    return _core.datas.btrans.get_sph_syms_no_cache()  # type: ignore


def get_sph_syms_prime_cache() -> list[str]:
    """获取球谐型基函数的对称性标签（预填充缓存）"""
    return _core.datas.btrans.get_sph_syms_prime_cache()  # type: ignore


def get_mix_syms() -> list[str]:
    """获取混合型基函数的对称性标签"""
    return _core.datas.btrans.get_mix_syms()  # type: ignore


def get_mix_syms_no_cache() -> list[str]:
    """获取混合型基函数的对称性标签（不使用缓存）"""
    return _core.datas.btrans.get_mix_syms_no_cache()  # type: ignore


def get_mix_syms_prime_cache() -> list[str]:
    """获取混合型基函数的对称性标签（预填充缓存）"""
    return _core.datas.btrans.get_mix_syms_prime_cache()  # type: ignore


def len2syms(n: int) -> list[str]:
    """根据基函数个数获取对应的对称性标签列表

    Args:
        n: 基函数个数
    """
    return _core.datas.btrans.len2syms(n)  # type: ignore


def len2syms_no_cache(n: int) -> list[str]:
    """根据基函数个数获取对称性标签（不使用缓存）"""
    return _core.datas.btrans.len2syms_no_cache(n)  # type: ignore


def len2syms_prime_cache(n: int) -> list[str]:
    """根据基函数个数获取对称性标签（预填充缓存）"""
    return _core.datas.btrans.len2syms_prime_cache(n)  # type: ignore


def lmn2sym(lmn: list[int]) -> str:
    """由量子数 (l, m, n) 获取对称性标签

    Args:
        lmn: 量子数三元组 [l, m, n]
    """
    return _core.datas.btrans.lmn2sym(lmn)  # type: ignore


def sym2ang(sym: str) -> int:
    """由对称性标签获取角动量量子数

    Args:
        sym: 对称性标签，如 "PX"
    """
    return _core.datas.btrans.sym2ang(sym)  # type: ignore


def sym2ang_no_cache(sym: str) -> int:
    """由对称性标签获取角动量量子数（不使用缓存）"""
    return _core.datas.btrans.sym2ang_no_cache(sym)  # type: ignore


def sym2ang_prime_cache(sym: str) -> int:
    """由对称性标签获取角动量量子数（预填充缓存）"""
    return _core.datas.btrans.sym2ang_prime_cache(sym)  # type: ignore


def sym2lmn(sym: str) -> list[int]:
    """由对称性标签获取量子数 (l, m, n)

    Args:
        sym: 对称性标签，如 "PX"
    """
    return _core.datas.btrans.sym2lmn(sym)  # type: ignore


def sym2lmn_no_cache(sym: str) -> list[int]:
    """由对称性标签获取量子数（不使用缓存）"""
    return _core.datas.btrans.sym2lmn_no_cache(sym)  # type: ignore


def sym2lmn_prime_cache(sym: str) -> list[int]:
    """由对称性标签获取量子数（预填充缓存）"""
    return _core.datas.btrans.sym2lmn_prime_cache(sym)  # type: ignore


def sym2stp(sym: str) -> int:
    """由对称性标签获取起始索引

    Args:
        sym: 对称性标签
    """
    return _core.datas.btrans.sym2stp(sym)  # type: ignore


def get_fmat(
    from_smat: np.ndarray,
    from_forms: list[str],
    into_forms: list[str],
    angs: list[int],
) -> np.ndarray:
    """按给定的基函数形式做基函数变换

    Args:
        from_smat: 原始矩阵
        from_forms: 原始矩阵每行对应的基函数形式
        into_forms: 目标矩阵每行对应的基函数形式
        angs: 每行对应的角动量量子数
    """
    return _core.datas.btrans.get_fmat(from_smat, from_forms, into_forms, angs)  # type: ignore