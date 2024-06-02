"""
为Fortran函数写python接口，方便调用
"""
from typing import Any
import ctypes as ct
from ctypes import c_long,c_float,POINTER
from pathlib import Path
import numpy as np
import os

def trans_dtype(paras:list):
    """转换数据类型"""
    targs=[]
    types=[]
    for para in paras:
        if type(para)==int:
            targs.append(c_long(para))
            types.append(c_long)
        elif type(para)==float:
            targs.append(c_float(para))
            types.append(c_float)
        elif type(para)==np.ndarray:
            # print(para.dtype,para.dtype in ['float32','float64'])
            if para.dtype in ['float32','float64']:
                if para.dtype!=np.float32:para=para.astype(np.float32) # 对于输出来说，不能进行类型转换，否则对原本对象的引用会改变
                targs.append(para.ctypes.data_as(POINTER(c_float)))
                types.append(POINTER(c_float))
            elif para.dtype in ['int32','int64']:
                if para.dtype!=np.int32:para=para.astype(np.int32) # 对于输出来说，不能进行类型转换，否则对原本对象的引用会改变
                targs.append(para.ctypes.data_as(POINTER(c_long)))
                types.append(POINTER(c_long))
            else:
                raise TypeError("Unsupported type")
        else:
            raise TypeError("Unsupported type")
    return targs,types

os.add_dll_directory(r"D:\program\mingw64\bin")

path = str(Path(__file__).parent / "flib.dll")
flib = ct.CDLL(path)

flib.add_.argtypes = [ct.c_float, ct.c_float, ct.POINTER(ct.c_float)]


def add_(x: float, y: float):
    res = ct.c_float()
    flib.add_(ct.c_float(x), ct.c_float(y), ct.byref(res))
    return res.value


# 抄的例子
flib.sum2_.argtypes = [ct.POINTER(ct.c_float)]
flib.sum2_.restype = ct.c_float


def sum2(a: float):
    a = ct.byref(ct.c_float(a))
    b = flib.sum2_(a)
    return b


# 抄的例子，矩阵运算
flib.double_array_.argtypes = [ct.POINTER(ct.c_float), ct.c_long]


def double_array():
    # Create a double array, pass it to Fotran as a pointer
    x = np.ones((3, 3), order="F")
    x_ptr = x.ctypes.data_as(ct.POINTER(ct.c_float))

    # Call function
    rint = flib.double_array_(x_ptr, ct.c_long(3))
    return x


flib.grid_pos_.argtypes = [ct.c_long, ct.c_long, ct.c_long, ct.POINTER(ct.c_float)]


def grid_pos(Nx: int, Ny: int, Nz: int):
    pos = np.zeros((Nx * Ny * Nz, 3), order="F")
    pos_ptr = pos.ctypes.data_as(ct.POINTER(ct.c_float))
    flib.grid_pos_(ct.c_long(Nx), ct.c_long(Ny), ct.c_long(Nz), pos_ptr)
    return pos


flib.same_array_.argtypes = [ct.c_long, ct.c_long, ct.POINTER(ct.c_float)]

def same_array():
    row, col = 2, 3
    pos = np.array([[1, 2, 3], [4, 5, 6]], dtype=float, order="F")

    pos_ptr = pos.ctypes.data_as(ct.POINTER(ct.c_float))
    flib.same_array_(ct.c_long(row), ct.c_long(col), pos_ptr)
    return pos


flib.gtf_.argtypes = [
    c_float, # alp
    c_long, # np
    ct.POINTER(ct.c_float), # pos
    ct.POINTER(ct.c_float), # r2
    ct.POINTER(ct.c_long), # lmn
    ct.POINTER(ct.c_float), # val
]

def gtf(exp: float, pos: np.ndarray, R2:np.ndarray, lmn: np.ndarray) -> np.ndarray:
    pos_ptr = pos.ctypes.data_as(ct.POINTER(ct.c_float))
    npos = pos.shape[0]
    R2_ptr = R2.ctypes.data_as(ct.POINTER(ct.c_float))
    lmn_ptr = lmn.ctypes.data_as(ct.POINTER(ct.c_long))
    wfn = np.zeros(npos)
    val_ptr = wfn.ctypes.data_as(ct.POINTER(ct.c_float))
    flib.gtf_(c_float(exp), c_long(npos), pos_ptr, R2_ptr, lmn_ptr, val_ptr)
    return wfn

def molDens(
        ngrid:int,
        grids:np.ndarray[(Any,3),float],
        nmat:int,
        cords:np.ndarray[(Any,3),float],
        nobt:int,
        CM:np.ndarray[float],
        ncgs:np.ndarray[int], # 每一个原子轨道的基组收缩数量
        cmax:int,
        oalps:np.ndarray[float],
        ocoes:np.ndarray[float],
        lmns:np.ndarray[int]
):
    """计算分子电子密度"""
    # print(grids[:3,:])
    paras=[ngrid,grids,nmat,cords,nobt,CM,ncgs,cmax,oalps,ocoes,lmns]
    iparas,itypes=trans_dtype(paras)

    dens=np.zeros(ngrid,dtype=np.float32)
    oparas,otypes=trans_dtype([dens])

    flib.moldens_.argtypes=itypes+otypes
    fparas=iparas+oparas
    # print(ncgs)
    # print(lmns)
    flib.moldens_(*fparas)
    return dens

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

    paras=[atm,nGrid,atmGrid,atmWeit,natm,atmPos,atmRad,atmDis]
    iparas,itypes=trans_dtype(paras)

    a2mGrid=np.zeros_like(atmGrid,dtype=np.float32)
    a2mWeit=np.zeros_like(atmWeit,dtype=np.float32)

    oparas,otypes=trans_dtype([a2mGrid,a2mWeit])

    if flib.a2mWeight_.argtypes is None:
        flib.a2mWeight_.argtypes=itypes+otypes

    flib.a2mWeight_(*(iparas+oparas))
    return a2mGrid,a2mWeit

import ctypes

# flib.test2_.argtypes = [ctypes.c_long] * 5