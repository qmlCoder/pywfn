"""
为Fortran函数写python接口，方便调用
"""
from typing import Any
import ctypes as ct
from ctypes import c_int,c_double,POINTER,byref
from pathlib import Path
import numpy as np
import os

ftype=np.float64

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
            # print(para.dtype,para.dtype in ['float32','float64'])
            if para.dtype in ['float32','float64']:
                if para.dtype!=ftype:para=para.astype(ftype) # 对于输出来说，不能进行类型转换，否则对原本对象的引用会改变
                targs.append(para.ctypes.data_as(POINTER(c_double)))
                types.append(POINTER(c_double))
            elif para.dtype in ['int32','int64']:
                if para.dtype!=np.int32:para=para.astype(np.int32) # 对于输出来说，不能进行类型转换，否则对原本对象的引用会改变
                targs.append(para.ctypes.data_as(POINTER(c_int)))
                types.append(POINTER(c_int))
            else:
                raise TypeError("Unsupported type")
        else:
            raise TypeError("Unsupported type")
    return targs,types

from pywfn import config
os.add_dll_directory(rf"{config.ROOT_LIBS}") # 添加动态链接库目录
flib = ct.CDLL(f'{config.ROOT_LIBS}/flib.dll')

flib.add_.argtypes = [ct.c_double, ct.c_double, ct.POINTER(ct.c_double)]


def add_(x: float, y: float):
    res = ct.c_double()
    flib.add_(ct.c_double(x), ct.c_double(y), ct.byref(res))
    return res.value


# 抄的例子
flib.sum2_.argtypes = [ct.POINTER(ct.c_double)]
flib.sum2_.restype = ct.c_double


def sum2(a: float):
    # a_ref = ct.byref(ct.c_double(a))
    a_ref = ct.c_double(a)
    b = flib.sum2_(a_ref)
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

def gtf(alp: float,ngrid:int, grids: np.ndarray,coord:np.ndarray,l:int,m:int,n:int) -> np.ndarray:
    assert len(coord.shape)==1,"坐标的维度为1"
    paras=[alp,ngrid,grids,coord,l,m,n]
    iparas,itypes=trans_dtype(paras)

    wfn=np.zeros(ngrid,dtype=ftype)
    oparas,otypes=trans_dtype([wfn])

    flib.gtf_.argtypes=itypes+otypes
    fparas=iparas+oparas
    flib.gtf_(*fparas)
    return wfn

def cgf(cmax:int,nc:int,alps:np.ndarray,coes:np.ndarray,
        ngrid:np.ndarray,grids:np.ndarray,coord:np.ndarray,l:int,m:int,n:int):
    assert len(coord.shape)==1,"坐标的维度为1"
    paras=[cmax,nc,alps,coes,ngrid,grids,coord,l,m,n]
    iparas,itypes=trans_dtype(paras)

    wfn=np.zeros(ngrid,dtype=ftype)
    oparas,otypes=trans_dtype([wfn])

    flib.cgf_.argtypes=itypes+otypes
    fparas=iparas+oparas
    flib.cgf_(*fparas)
    return wfn

Array=np.ndarray

def cgfs(ngrid:int,grids:Array,nmat:int,cords:Array,cmax:int,ncgs:Array,alpl:Array,coel:Array,lmns:Array):
    """
    计算所有原子轨道波函数
    """
    paras=[ngrid,grids,nmat,cords,cmax,ncgs,alpl,coel,lmns]
    iparas,itypes=trans_dtype(paras)

    wfns=np.zeros(shape=(nmat,ngrid),dtype=ftype)
    oparas,otypes=trans_dtype([wfns])

    flib.cgfs_.argtypes=itypes+otypes
    fparas=iparas+oparas
    flib.cgfs_(*fparas)

    return wfns

def molDens(ngrid:int,nmat:int,nobt:int,CM:np.ndarray,wfns:np.ndarray):
    """计算分子电子密度"""
    # print(grids[:3,:])
    paras=[ngrid,nmat,nobt,CM,wfns]
    iparas,itypes=trans_dtype(paras)

    dens=np.zeros(ngrid,dtype=ftype)
    oparas,otypes=trans_dtype([dens])

    flib.moldens_.argtypes=itypes+otypes
    fparas=iparas+oparas
    flib.moldens_(*fparas)
    return dens

# def atmDens()

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

    a2mGrid=np.zeros_like(atmGrid,dtype=ftype)
    a2mWeit=np.zeros_like(atmWeit,dtype=ftype)
    total=c_int()
    oparas,otypes=trans_dtype([a2mGrid,a2mWeit])

    # if flib.a2mWeight_.argtypes is None:
    flib.a2mWeight_.argtypes=itypes+otypes+[POINTER(c_int)]
    flib.a2mWeight_(*(iparas+oparas+[byref(total)]))
    # print(f'{total=},{len(atmWeit)}')
    return a2mGrid[:total.value,:],a2mWeit[:total.value]

def eleMat(nmat:int,nobt:int,CM:np.ndarray,SM:np.ndarray)->np.ndarray:
    """计算电子数量矩阵"""
    assert CM.shape==(nmat,nobt),'CM.shape must be (nmat,nobt)'
    assert SM.shape==(nmat,nmat),'SM.shape must be (nmat,nmat)'
    paras=[nmat,nobt,CM,SM]
    iparas,itypes=trans_dtype(paras)
    
    NM=np.zeros((nmat,nobt),dtype=ftype)
    oparas,otypes=trans_dtype([NM])
    flib.eleMat_.argtypes=itypes+otypes
    fparas=iparas+oparas
    flib.eleMat_(*fparas)
    return NM