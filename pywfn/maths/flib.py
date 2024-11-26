"""
为Fortran函数写python接口，方便调用
Fortran函数的流程
接受输入参数 ipts
转换输入参数 iparas,itypes
定义输出参数 outs
转换输出参数 oparas,otypes
设置参数类型 fargs,argtypes

参数传递顺序：
输入在前，输出在后
整数在前，小数在后
数值在前，数组在后
"""
from typing import Any
import ctypes as ct
from ctypes import c_int,c_double,POINTER,byref
from pathlib import Path
import numpy as np
import os
from pywfn.utils import chkArray

ftype=np.float64 # 使用的浮点数类型
itype=np.int32   # 使用的整数类型
Array=np.ndarray

def trans_dtype(paras:list):
    fparas=[]
    ftypes=[]
    for para in paras:
        if type(para)==int:
            fparas.append(c_int(para))
            ftypes.append(c_int)
        elif type(para)==float:
            fparas.append(c_double(para))
            ftypes.append(c_double)
        elif type(para)==np.ndarray:
            if para.dtype in ['int32','int64']:
                if para.dtype!=itype:
                    para=para.astype(itype) # 对于输出来说，不能进行类型转换，否则对原本对象的引用会改变
                # print(para.dtype)
                fparas.append(para.ctypes.data_as(POINTER(c_int)))
                ftypes.append(POINTER(c_int))
            elif para.dtype in ['float32','float64']:
                if para.dtype!=ftype:para=para.astype(ftype) # 对于输出来说，不能进行类型转换，否则对原本对象的引用会改变
                fparas.append(para.ctypes.data_as(POINTER(c_double)))
                ftypes.append(POINTER(c_double))
            else:
                raise TypeError("Unsupported type")
        else:
            raise TypeError(f"Unsupported type: {type(para)}")
    return fparas,ftypes

def call_flib(func:str,ipts:list,outs:list):
    """调用Fortran函数

    Args:
        func (str): 函数名
        ipts (list): 输入参数
        outs (list): 输出参数
    """
    # print(func)
    iparas,itypes=trans_dtype(ipts)
    oparas,otypes=trans_dtype(outs)
    ftypes=itypes+otypes
    flib[func].argtypes=ftypes # 指定函数参数的类型
    fparas=iparas+oparas  # 指定函数参数的值
    flib[func](*fparas)

from pywfn import config
print('动态链接库目录',config.ROOT_LIBS)
if os.name=='nt': # Windows系统
    print('当前系统:windows')
    os.add_dll_directory(rf"{config.ROOT_LIBS}") # 添加动态链接库目录
    flib = ct.CDLL(f'{config.ROOT_LIBS}/flib.dll')
elif os.name=='posix': # Linux系统
    print('当前系统:linux')
    os.environ['LD_LIBRARY_PATH'] = rf"{config.ROOT_LIBS}" + os.environ.get('LD_LIBRARY_PATH', '')
    flib = ct.CDLL(f'{config.ROOT_LIBS}/flib.so')
else:
    raise OSError("Unsupported OS")
call_flib('info_',[],[])

def grid_pos(Nx: int, Ny: int, Nz: int)->np.ndarray:
    """生成格点数据

    Args:
        Nx (int): x方向数量
        Ny (int): y方向数量
        Nz (int): z方向数量

    Returns:
        np.ndarray: 格点数据
    """ 
    ipts=[Nx,Ny,Nz]
    pos = np.zeros((Nx * Ny * Nz, 3))
    outs=[pos]

    call_flib('grid_pos_',ipts,outs)
    return pos

def gtf(alp: float,ngrid:int, grids: np.ndarray,coord:np.ndarray,l:int,m:int,n:int) -> np.ndarray:
    assert chkArray(coord,[3]),"形状或类型不匹配"
    ipts=[alp,ngrid,grids,coord,l,m,n]
    wfn=np.zeros(ngrid,dtype=ftype)
    outs=[wfn]
    call_flib('gtf_',ipts,outs)
    return wfn


def cgf(cmax:int,nc:int,alps:np.ndarray,coes:np.ndarray,
        ngrid:np.ndarray,grids:np.ndarray,coord:np.ndarray,l:int,m:int,n:int):
    assert len(coord.shape)==1,"坐标的维度为1"
    wfn=np.zeros(ngrid,dtype=ftype)

    ipts=[cmax,nc,alps,coes,ngrid,grids,coord,l,m,n]
    outs=[wfn]
    call_flib('cgf_',ipts,outs)
    return wfn

def atoWfns(
        ngrid:int,
        grids:Array,
        nmat:int,
        cords:Array,
        cmax:int,
        ncgs:Array,
        alpl:Array,
        coel:Array,
        lmns:Array):
    """
    计算所有原子轨道波函数
    """
    assert chkArray(grids,[ngrid,]),"形状不匹配"
    wfns=np.zeros(shape=(nmat,ngrid),dtype=ftype)
    ipts=[ngrid,grids,nmat,cords,cmax,ncgs,alpl,coel,lmns]
    outs=[wfns]
    call_flib('atoWfns_',ipts,outs)
    return wfns

def molDens(ngrid:int,nmat:int,nobt:int,CM:np.ndarray,wfns:np.ndarray):
    """计算分子电子密度"""
    # print(grids[:3,:])
    paras=[ngrid,nmat,nobt,CM,wfns]
    dens=np.zeros(ngrid,dtype=ftype)
    call_flib('moldens_',paras,[dens])
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
    """计算dft格点权重

    Args:
        atm (int): 原子索引
        nGrid (int): 格点数量
        atmGrid (np.ndarray): 原子dft格点
        natm (int): 原子数量
        atmPos (np.ndarray): 原子坐标
        atmRad (np.ndarray): 原子半径
        atmDis (np.ndarray): 原子距离

    Returns:
        _type_: _description_
    """
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

def nucPotential(cords:Array,nucs:Array,xyzs:Array):
    """计算势能"""
    assert chkArray(cords,[None,3]),"形状不匹配"
    assert chkArray(nucs,[None,]),"形状不匹配"
    assert chkArray(xyzs,[None,3]),"形状不匹配"
    ncord=cords.shape[0]
    natm=len(nucs)
    paras=[ncord,cords,natm,nucs,xyzs]
    
    vals=np.zeros(ncord,dtype=ftype)
    call_flib('nucPotential_',paras,[vals])

    return vals

def elePotential(cords:Array,grids:Array,weits:Array,dens:Array):
    """计算势能"""
    assert chkArray(cords,[None,3]),"形状不匹配"
    assert chkArray(grids,[None,3]),"形状不匹配"
    assert chkArray(weits,[None,]),"形状不匹配"
    assert chkArray(dens,[None,]),"形状不匹配"
    ncord=cords.shape[0]
    ngrid=grids.shape[0]
    paras=[ncord,cords,ngrid,grids,weits,dens]
    vals=np.zeros(ncord,dtype=ftype)
    call_flib('elePotential_',paras,[vals])
    return vals
