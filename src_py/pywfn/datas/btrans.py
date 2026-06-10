"""
记录基函数变换矩阵
http://sobereva.com/97 感谢卢天大大
"""

import numpy as np

from pywfn import _core

DMAT: np.ndarray = _core.datas.btrans.DMAT()  # type: ignore
FMAT: np.ndarray = _core.datas.btrans.FMAT()  # type: ignore
GMAT: np.ndarray = _core.datas.btrans.GMAT()  # type: ignore
HMAT: np.ndarray = _core.datas.btrans.HMAT()  # type: ignore

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
