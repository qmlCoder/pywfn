"""
公用的函数可以移到此处，以较少代码量
"""
from pywfn.base import Atom,Mole
import numpy as np
from pywfn import maths


def atomIdxs(atoms:list[Atom])->list[int]:
    return [atom.idx for atom in atoms]


def get_ects(mol:Mole,obts:list[int],CM:np.ndarray)->list[int]:
    """
    计算电子数量,如果不指定则计算所有的π电子
    如果指定方向，则计算方向电子
    每个原子都有应该根据自己的法向量方向来计算(当法向量都相同时相当于都一样)
    """
    # 计算密度矩阵
    PM=maths.CM2PM(CM,obts,mol.oE) # 计算密度矩阵
    # PS=PM*mol.SM
    # PSS=PS.sum(axis=0)
    # 矩阵乘法的迹的加和=矩阵对应元素乘积之和
    PS=PM@mol.SM
    EV=np.diagonal(PS)
    elects=[]
    for atom in mol.atoms:
        a1,a2=atom.obtBorder
        elect=EV[a1:a2].sum()
        elects.append(elect)
    return elects

def atomValueStr(mol:Mole,satoms:list[int],values:list[float]):
    resStr=''
    if values is None:return '非开窍层无法计算自旋!!'
    
    for i,atom in enumerate(mol.atoms):
        if atom.idx in satoms:
            idx=satoms.index(atom.idx)
            value=values[idx]
        else:
            continue
        resStr+=f'{atom.idx:<2}{atom.symbol:>2}{value:>15.8f}\n'

    return resStr


def mulliken(
        PM:np.ndarray,
        SM:np.ndarray,
        atmuls:list[tuple[int,int]], # 每个原子在基函数中的上下限
        )->np.ndarray:
    """
    计算目录mulliken电荷
    num：是否只保留电子数
    """
    # 矩阵乘法的迹的加和=矩阵对应元素乘积之和
    PS=PM@SM
    EV=np.diagonal(PS) # 矩阵的对角元素
    elects=np.zeros(len(atmuls))
    for a,(u,l) in enumerate(atmuls):
        elect=EV[u:l].sum()
        elects[a]=elect
    return elects

def lowdin(
        PM:np.ndarray,
        SM:np.ndarray,
        atmuls:list[tuple[int,int]], # 每个原子在基函数中的上下限
    )->np.ndarray:
    """
    计算每个原子的lowdin电荷
    """
    # 计算矩阵的1/2次方
    v, Q=np.linalg.eig(SM)
    V=np.diag(v)
    V_=V**0.5
    Q_=np.linalg.inv(Q)
    SM_half=Q@(V_@Q_)
    SPS=SM_half@(PM@SM_half)
    EV=np.diagonal(SPS)
    elects=np.zeros(len(atmuls))
    for a,(u,l) in enumerate(atmuls):
        elect=EV[u:l].sum()
        elects[a]=elect
    return elects