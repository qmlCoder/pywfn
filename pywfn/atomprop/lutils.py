"""
公用的函数可以移到此处，以较少代码量
"""
from pywfn.base import Atom,Mol
import numpy as np
from pywfn import maths


def atomIdxs(atoms:list[Atom])->list[int]:
    return [atom.idx for atom in atoms]

def get_vects(mol:Mol,atoms:list[int]=None):
    """获取指定原子的法向量或与相邻的原子的法向量"""
    atoms:list[Atom]=[mol.atom(a) for a in atoms]
    vects=[]
    for atom in atoms:
        vect=atom.get_Normal()
        if vect is None:
            idxn=maths.search_sp2(atom.idx,mol)
            vect=mol.atom(idxn).get_Normal()
        vects.append(vect)
    return vects

def get_ects(mol:Mol,obts,CM_)->list[int]:
    """
    计算电子数量,如果不指定则计算所有的π电子
    如果指定方向，则计算方向电子
    每个原子都有应该根据自己的法向量方向来计算(当法向量都相同时相当于都一样)
    """
    # 计算密度矩阵
    PM_=maths.CM2PM(CM_,obts,mol.oE) # 计算密度矩阵
    PS=PM_*mol.SM
    PSS=PS.sum(axis=0)
    elects=[]
    for atom in mol.atoms:
        a_1,a_2=atom.obtBorder
        electron=np.sum(PSS[a_1:a_2])
        elects.append(electron)
    return elects

def atomValueStr(mol:Mol,satoms:list[int],values:list[float]):
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