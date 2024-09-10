from pywfn.base import Mol
from pywfn.maths import vector_angle

import numpy as np
import re

def dihedralAngle(mol:Mol,idxs:list[int]):
    """计算二面角"""
    a,b,c,d=[mol.atom(idx) for idx in idxs]
    vba=a.coord-b.coord
    vbc=c.coord-b.coord
    vcd=d.coord-c.coord
    vcb=b.coord-c.coord
    vi=np.cross(vba,vbc)
    vj=np.cross(vcb,vcd)
    angle=vector_angle(vi,vj)
    return angle

def projCM(mol:Mol,obts:list[int],atms:list[int],dirs:list[np.ndarray],
           akeep:bool,lkeep:bool,akeeps=None,keeps:str|None=None)->np.ndarray:
    """
    获取投影后的系数矩阵
    atms:需要投影的原子
    obts:需要投影的轨道
    dirs:投影到的方向 atms和dirs的长度必须相同
    akeep:其它原子系数是否保留 keep other atoms
    lkeep:其它价层系数是否保留 keep other layer
    akeeps:额外保留的原子，不进行投影但是保留
    lkeeps:额外保留的轨道，可以使用正则表达式匹配
    """
    assert isinstance(dirs,list),"方向向量需为列表"
    assert len(atms)==len(dirs),"原子和方向数量不同"
    if akeep:
        CMp=np.copy(mol.CM)
    else:
        CMp=np.zeros_like(mol.CM,dtype=float) #新的系数矩阵
        
    for a,(atom,vect) in enumerate(zip(atms,dirs)):
        atom=mol.atom(atom)
        nebNum=len(atom.neighbors)
        u,l=atom.obtBorder
        syms=mol.obtSyms[u:l]
        
        pIdx=[i for i,s in enumerate(syms) if 'P' in s]
        if len(pIdx)==0:continue # 没有p轨道则跳过
        for o,obt in enumerate(obts):
            if lkeep:
                Co=mol.CM[u:l,obt].copy()
            else:
                Co=np.zeros(len(syms)) #根据是否P轨道之外的保留还是0由不同的选择
            if keeps is not None: # 其它保留的层
                idxs=[i for i,s in enumerate(syms) if re.match(keeps,s)]
                Cos=atom.obtCoeffs[idxs,obt].copy()
                Co[idxs]=Cos
            Cop=atom.get_pProj(vect,obt)
            Co[pIdx]=np.concatenate(Cop)
            CMp[u:l,obt]=Co.copy()
    if akeeps is not None: # 保留一些指定的原子的系数
        for a in akeeps:
            atom=mol.atom(a)
            u,l=atom.obtBorder
            CMp[u:l]=mol.CM[u:l].copy()
    return CMp