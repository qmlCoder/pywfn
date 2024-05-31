"""
为Fortran函数写python接口，方便调用
"""

import ctypes as ct
from ctypes import c_int, c_double,POINTER
from pathlib import Path
import numpy as np
import os

def trans_dtype(paras:list):
    """转换数据类型"""
    targs=[]
    types=[]
    for para in paras:
        if type(para)==int:
            targs.append(c_int(para))
            types.append(c_int)
        elif type(para)==float:
            targs.append(c_double(para))
            types.append(c_double)
        elif type(para)==np.ndarray:
            targs.append(para.ctypes.data_as(POINTER(c_double)))
            types.append(POINTER(c_double))
        else:
            raise TypeError("Unsupported type")
    return targs,types

os.add_dll_directory(r"D:\program\mingw64\bin")

path = str(Path(__file__).parent / "flib.dll")
flib = ct.CDLL(path)

flib.add_.argtypes = [ct.c_double, ct.c_double, ct.POINTER(ct.c_double)]


def add_(x: float, y: float):
    res = ct.c_double()
    flib.add_(ct.c_double(x), ct.c_double(y), ct.byref(res))
    return res.value


# 抄的例子
flib.sum2_.argtypes = [ct.POINTER(ct.c_double)]
flib.sum2_.restype = ct.c_double


def sum2(a: float):
    a = ct.byref(ct.c_double(a))
    b = flib.sum2_(a)
    return b


# 抄的例子，矩阵运算
flib.double_array_.argtypes = [ct.POINTER(ct.c_double), ct.c_int]


def double_array():
    # Create a double array, pass it to Fotran as a pointer
    x = np.ones((3, 3), order="F")
    x_ptr = x.ctypes.data_as(ct.POINTER(ct.c_double))

    # Call function
    rint = flib.double_array_(x_ptr, ct.c_int(3))
    return x


flib.grid_pos_.argtypes = [ct.c_int, ct.c_int, ct.c_int, ct.POINTER(ct.c_double)]


def grid_pos(Nx: int, Ny: int, Nz: int):
    pos = np.zeros((Nx * Ny * Nz, 3), order="F")
    pos_ptr = pos.ctypes.data_as(ct.POINTER(ct.c_double))
    flib.grid_pos_(ct.c_int(Nx), ct.c_int(Ny), ct.c_int(Nz), pos_ptr)
    return pos


flib.same_array_.argtypes = [ct.c_int, ct.c_int, ct.POINTER(ct.c_double)]

def same_array():
    row, col = 2, 3
    pos = np.array([[1, 2, 3], [4, 5, 6]], dtype=float, order="F")

    pos_ptr = pos.ctypes.data_as(ct.POINTER(ct.c_double))
    flib.same_array_(ct.c_int(row), ct.c_int(col), pos_ptr)
    return pos


flib.gtf_.argtypes = [
    c_double, # alp
    c_int, # np
    ct.POINTER(ct.c_double), # pos
    ct.POINTER(ct.c_double), # r2
    ct.POINTER(ct.c_int), # lmn
    ct.POINTER(ct.c_double), # val
]

def gtf(exp: float, pos: np.ndarray, R2:np.ndarray, lmn: np.ndarray) -> np.ndarray:
    pos_ptr = pos.ctypes.data_as(ct.POINTER(ct.c_double))
    npos = pos.shape[0]
    R2_ptr = R2.ctypes.data_as(ct.POINTER(ct.c_double))
    lmn_ptr = lmn.ctypes.data_as(ct.POINTER(ct.c_int))
    wfn = np.zeros(npos)
    val_ptr = wfn.ctypes.data_as(ct.POINTER(ct.c_double))
    print(R2.shape)
    flib.gtf_(c_double(exp), c_int(npos), pos_ptr, R2_ptr, lmn_ptr, val_ptr)
    return wfn



flib.a2mWeight_.argtypes=[
    ct.c_int, # atm
    ct.c_int, # nGrid
    ct.POINTER(c_double), # atmGrid
    ct.POINTER(c_double), # atmWeit
    ct.c_int, # natm
    ct.POINTER(c_double), # atmPos
    ct.POINTER(c_double), # atmRad
    ct.POINTER(c_double), # atmDis
    ct.POINTER(c_double), # a2mGrid
    ct.POINTER(c_double), # a2mWeit
]
def a2mWeight(
        atm:int,
        nGrid:int,
        atmGrid:np.ndarray, # 原子的网格点坐标
        atmWeit:np.ndarray,
        natm:int,
        atmPos:np.ndarray,
        atmRad:np.ndarray,
        atmDis:np.ndarray,
        ):
    """
    """
    a2mGrid=np.zeros_like(atmGrid)
    a2mWeit=np.zeros_like(atmWeit)
    flib.a2mWeight_(
        ct.c_int(atm),
        ct.c_int(nGrid),
        atmGrid.ctypes.data_as(ct.POINTER(c_double)),
        atmWeit.ctypes.data_as(ct.POINTER(c_double)),
        ct.c_int(natm),
        atmPos.ctypes.data_as(ct.POINTER(c_double)),
        atmRad.ctypes.data_as(ct.POINTER(c_double)),
        atmDis.ctypes.data_as(ct.POINTER(c_double)),
        a2mGrid.ctypes.data_as(ct.POINTER(c_double)),
        a2mWeit.ctypes.data_as(ct.POINTER(c_double)),
    )
    return a2mGrid,a2mWeit
    # flib.a2mWeight_(atm,nGrid,atmGrid,atmWeit,natm,atmPos,atmRad,atmDis,a2mGrid,a2mWeit)

