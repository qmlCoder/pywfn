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
import ctypes
import ctypes as ct
from ctypes import c_int,c_double,POINTER,byref
# import faulthandler

from pathlib import Path
import numpy as np
import os
from pywfn.utils import chkArray

ftype=np.float64 # 使用的浮点数类型
itype=np.int32   # 使用的整数类型
Array=np.ndarray


# faulthandler.enable()

def trans_dtype(paras:list,intent='in'): # 将参数列表转为ctypes需要的类型
    fparas=[]
    ftypes=[]
    for para in paras:
        if isinstance(para, int): # 如果是整数输入，则转换为c_int
            if intent=='in':
                fparas.append(byref(c_int(para)))
            ftypes.append(POINTER(c_int))
            # print('整数标量',para)
        elif isinstance(para, float):
            if intent=='in':
                fparas.append(byref(c_double(para)))
            ftypes.append(POINTER(c_double))
            # print('小数标量',para)
        elif isinstance(para, np.ndarray):
            # if intent=='in':
            # para=np.copy(para) # 防止对原始数据的修改
            if para.dtype in ['int32','int64']:
                if para.dtype!=itype:
                    para=para.astype(itype) # 对于输出来说，不能进行类型转换，否则对原本对象的引用会改变
                fparas.append(para.ctypes.data_as(POINTER(c_int)))
                ftypes.append(POINTER(c_int))
                # print('整数数组',para.shape)
            elif para.dtype in ['float32','float64']:
                if para.dtype!=ftype:
                    para=para.astype(ftype) # 对于输出来说，不能进行类型转换，否则对原本对象的引用会改变
                fparas.append(para.ctypes.data_as(POINTER(c_double)))
                ftypes.append(POINTER(c_double))
                # print('小数数组',para.shape)
            else:
                raise TypeError("Unsupported type")
        elif str(type(para))=="<class 'CArgObject'>":
            fparas.append(para)
            ftypes.append(POINTER(c_int))
        else:
            # raise TypeError(f"Unsupported type: {type(para)}")
            print(f"类型转换失败: {type(para)}")
    return fparas,ftypes

def call_flib(func:str,ipts:list,outs:list):
    """调用Fortran函数

    Args:
        func (str): 函数名
        ipts (list): 输入参数
        outs (list): 输出参数
    """
    # print(func)
    iparas,itypes=trans_dtype(ipts,'in')
    oparas,otypes=trans_dtype(outs,'out')
    ftypes=itypes+otypes
    flib[func].argtypes=ftypes # 指定函数参数的类型
    fparas=iparas+oparas  # 指定函数参数的值
    flib[func](*fparas)

from pywfn import config
# print('动态链接库目录',config.ROOT_LIBS)
if os.name=='nt': # Windows系统
    # print(f'当前系统:windows, 动态链接库目录{config.ROOT_LIBS}')
    os.add_dll_directory(rf"{config.ROOT_LIBS}") # 添加动态链接库目录
    flib = ct.CDLL(f'{config.ROOT_LIBS}/flib.dll')
    flib['init_']()
elif os.name=='posix': # Linux系统
    # print('当前系统:linux') 你需要使用gfortran自行编译动态链接库，
    # 使用conda安装gcc gxx gfortran: conda install -c conda-forge gcc=14 gxx=14 gfortran=14
    # 编译命令为 gfortran -shared -ffree-line-length-none -fopenmp -fPIC data.f90 march.f90 flib.f90 -o flib.so
    os.environ['LD_LIBRARY_PATH'] = rf"{config.ROOT_LIBS}" + os.environ.get('LD_LIBRARY_PATH', '')
    flib = ct.CDLL(f'{config.ROOT_LIBS}/flib.so')
    flib.init_()
else:
    raise OSError("Unsupported OS")
# call_flib('info_',[],[])

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

def gtf(alp: float,ngrid:int, grids: np.ndarray,coord:np.ndarray,l:int,m:int,n:int,level:int):
    assert chkArray(coord,[3,]),"形状或类型不匹配"
    ipts=[alp,ngrid,grids,coord,l,m,n,level]
    wfn0=np.zeros(shape=(ngrid,),dtype=ftype)
    wfn1=np.zeros(shape=(ngrid,3),dtype=ftype)
    wfn2=np.zeros(shape=(ngrid,3,3),dtype=ftype)
    outs=[wfn0,wfn1,wfn2]
    call_flib('gtf_',ipts,outs)
    return wfn0,wfn1,wfn2


def cgf(cmax:int,nc:int,alps:np.ndarray,coes:np.ndarray,
        ngrid:np.ndarray,grids:np.ndarray,coord:np.ndarray,l:int,m:int,n:int,level:int):
    assert len(coord.shape)==1,"坐标的维度为1"
    wfn0=np.zeros(shape=(ngrid,),dtype=ftype)
    wfn1=np.zeros(shape=(ngrid,3),dtype=ftype)
    wfn2=np.zeros(shape=(ngrid,3,3),dtype=ftype)

    ipts=[cmax,nc,alps,coes,ngrid,grids,coord,l,m,n,level]
    outs=[wfn0,wfn1,wfn2]
    call_flib('cgf_',ipts,outs)
    return wfn0,wfn1,wfn2

def atoWfns(
        ngrid:int,
        grids:Array,
        nmat:int,
        coords:Array, # 原子坐标
        cmax:int,
        ncgs:Array,
        alpl:Array,
        coel:Array,
        lmns:Array,
        level:int):
    """
    计算所有原子轨道波函数
    """
    assert chkArray(grids,[ngrid,3]),"形状不匹配"
    assert chkArray(coords,[nmat,3]),"形状不匹配"
    wfns0=np.zeros(shape=(nmat,ngrid),dtype=ftype)
    wfns1=np.zeros(shape=(nmat,ngrid,3))
    wfns2=np.zeros(shape=(nmat,ngrid,3,3))
    ipts=[ngrid,grids,nmat,coords,cmax,ncgs,alpl,coel,lmns,level]
    outs=[wfns0,wfns1,wfns2]
    call_flib('atoWfns_',ipts,outs)
    return wfns0,wfns1,wfns2

def obtDens(ngrid:int,nmat:int,nobt:int,CM:np.ndarray,wfns0:np.ndarray,wfns1:np.ndarray,wfns2:np.ndarray,level:int):
    paras=[ngrid,nmat,nobt,CM,wfns0,wfns1,wfns2,level]
    dens0=np.zeros(shape=(nobt,ngrid,),dtype=ftype)
    dens1=np.zeros(shape=(nobt,ngrid,3),dtype=ftype)
    dens2=np.zeros(shape=(nobt,ngrid,3,3),dtype=ftype)
    # print(dens0.size+dens1.size+dens2.size)
    call_flib('obtDens_',paras,[dens0,dens1,dens2])
    return dens0,dens1,dens2

def atoDens(ngrid:int,nmat:int,PM:np.ndarray,wfns0:np.ndarray,wfns1:np.ndarray,wfns2:np.ndarray,level:int):
    paras=[ngrid,nmat,PM,wfns0,wfns1,wfns2,level]
    dens0=np.zeros(shape=(nmat,ngrid,),dtype=ftype)
    dens1=np.zeros(shape=(nmat,ngrid,3),dtype=ftype)
    dens2=np.zeros(shape=(nmat,ngrid,3,3),dtype=ftype)
    # print(dens0.size+dens1.size+dens2.size)
    call_flib('atoDens_',paras,[dens0,dens1,dens2])
    return dens0,dens1,dens2

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
    # iparas,itypes=trans_dtype(paras)

    a2mGrid=np.zeros_like(atmGrid,dtype=ftype)
    a2mWeit=np.zeros_like(atmWeit,dtype=ftype)
    total=c_int()
    # oparas,otypes=trans_dtype([a2mGrid,a2mWeit])

    # # if flib.a2mWeight_.argtypes is None:
    # flib.a2mWeight_.argtypes=itypes+otypes+[POINTER(c_int)]
    # flib.a2mWeight_(*(iparas+oparas+[byref(total)]))
    # print(f'{total=},{len(atmWeit)}')
    call_flib('a2mWeight_',paras,[a2mGrid,a2mWeit,byref(total)])
    return a2mGrid[:total.value,:].copy(),a2mWeit[:total.value].copy()

def eleMat(nmat:int,nobt:int,CM:np.ndarray,SM:np.ndarray)->np.ndarray:
    """计算电子数量矩阵"""
    assert CM.shape==(nmat,nobt),'CM.shape must be (nmat,nobt)'
    assert SM.shape==(nmat,nmat),'SM.shape must be (nmat,nmat)'
    paras=[nmat,nobt,CM,SM]
    # iparas,itypes=trans_dtype(paras)
    
    NM=np.zeros((nmat,nobt),dtype=ftype)
    # oparas,otypes=trans_dtype([NM])
    # flib.eleMat_.argtypes=itypes+otypes
    # fparas=iparas+oparas
    # flib.eleMat_(*fparas)
    call_flib('eleMat_',paras,[NM])
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
    ncord=cords.shape[0]
    ngrid=grids.shape[0]
    paras=[ncord,cords,ngrid,grids,weits,dens]
    vals=np.zeros(ncord,dtype=ftype)
    call_flib('elePotential_',paras,[vals])
    return vals

def vertsShift(verts:Array)->Array:
    """顶点移动"""
    assert chkArray(verts,[None,3]),"形状不匹配"
    nvert=verts.shape[0]
    paras=[nvert]
    call_flib('vertsShift_',paras,[verts])
    return verts

def vertsMerge(old_verts:Array,thval:float):
    """顶点合并""" 
    assert chkArray(old_verts,[None,3]),f"形状不匹配,{old_verts.shape}"
    nvert=old_verts.shape[0]
    new_verts=np.zeros((nvert,3),dtype=ftype)
    faces=np.zeros(nvert,dtype=itype)
    vCount=c_int(0)
    fCount=c_int(0)
    call_flib('vertsMerge_',[thval,nvert,old_verts],[new_verts,faces,byref(vCount),byref(fCount)])
    # print(vCount.value)
    new_verts=new_verts[:vCount.value,:].copy()
    faces=faces[:fCount.value*3].copy()
    return new_verts,faces-1

def marchCube(grids:Array,values:Array,shape:list[int],isov:float):
    """顶点提取"""
    nx,ny,nz=shape
    voxel=np.zeros(shape=(nx,ny,nz,4),dtype=ftype) # 体素
    call_flib('grids2voxel_',[nx,ny,nz,grids,values],[voxel])
    # print(voxel)
    verts=np.zeros(shape=(nx*ny*nz*24,3),dtype=ftype) # 顶点
    count=c_int(0)
    call_flib('voxel2verts_',[nx,ny,nz,voxel,isov],[verts,byref(count)])
    return voxel,verts[:count.value]

def matInteg(atos:Array,coes:Array,alps:Array,lmns:Array,xyzs:Array):
    """计算重叠矩阵

    Args:
        atos (Array): 基函数对应的原子轨道
        coes (Array): 基函数的系数
        alps (Array): 基函数的指数
        lmns (Array): 基函数的角动量
        xyzs (Array): 基函数的坐标

    Returns:
        Array: 重叠矩阵
    """
    nato=len(xyzs)
    nbas=len(coes)
    assert chkArray(xyzs,[nato,3]),"形状不匹配"
    assert chkArray(lmns,[nbas,3]),"形状不匹配"
    # paras=[nato,nbas,atos,coes,alps,lmns,xyzs]
    # 将系数转为归一化之后的系数
    facs = [1., 1., 3.]
    for i in range(nbas):
        l,m,n=lmns[i]
        fac = facs[l]*facs[m]*facs[n]
        ang=l+m+n
        alp=alps[i]
        Nm=(2.*alp/np.pi)**0.75*np.sqrt((4.*alp)**ang/fac)
        coes[i]=Nm*coes[i]
        # print(f'{alp:>10.4f}->{Nm:>10.4f}|{coes[i]:>10.4f}')
    # print('atos',atos)
    paras=[nato,nbas,atos,coes,alps,lmns,xyzs]
    SM=np.zeros(shape=(nato,nato),dtype=ftype)
    # print('atos',atos.shape,atos)
    call_flib('matInteg_',paras,[SM])
    return SM

def lagIntpol(xs:Array,ys:Array,ts:Array):
    """拉格朗日插值"""
    nx=len(xs)
    nt=len(ts)
    vs=np.zeros_like(ts)
    call_flib('lagIntpol_',[nx,xs,ys,nt,ts],[vs])
    return vs