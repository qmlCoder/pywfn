from pywfn.base import Mol
from pywfn.maths import vector_angle

from pywfn.utils import chkArray

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

def projCM(mol:Mol,
           obts:list[int],
           atms:list[int],
           dirs:np.ndarray,
           akeep:bool,
           lkeep:bool,
           akeeps:list[int]|None=None,
           keeps:str|None=None)->np.ndarray:
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
    assert chkArray(dirs,[None,3]),"方向数组形状不对"
    assert len(atms)==len(dirs),"原子和方向数量不同"
    if akeep:
        CMp=np.copy(mol.CM)
    else:
        CMp=np.zeros_like(mol.CM,dtype=float) #新的系数矩阵
        
    for a,(atom,vect) in enumerate(zip(atms,dirs)):
        atom=mol.atom(atom)
        nebNum=len(atom.neighbors)
        u,l=atom.obtBorder
        syms=mol.atoSyms[u:l] #该原子的轨道符号
        
        pIdx=[i for i,s in enumerate(syms) if 'P' in s] # p轨道的索引
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
            Cop=np.concatenate(Cop)
            Co[pIdx]=Cop
            CMp[u:l,obt]=Co.copy()
    if akeeps is not None: # 保留一些指定的原子的系数
        for a in akeeps:
            atom=mol.atom(a)
            u,l=atom.obtBorder
            CMp[u:l]=mol.CM[u:l].copy()
    return CMp

def hmo(mol:Mol)->tuple[np.ndarray,np.ndarray,np.ndarray,np.ndarray]:
    """求解休克尔分子轨道法

    Args:
        mol (Mol): 需要求解的分子

    Returns:
        tuple[np.ndarray,np.ndarray,np.ndarray]: 距离矩阵，能量，系数矩阵
    """
    from pywfn.atomprop import direction
    dirCaler=direction.Calculator(mol)
    atomBases=dirCaler.hmoBases()
    # 1.建立键连矩阵
    atms=mol.heavyAtoms # 重原子列表
    natm=len(atms) # 重原子数量
    BM=np.zeros(shape=(natm,natm)) # 键连矩阵
    # DM=np.zeros(shape=(natm,natm)) # 距离矩阵
    for i,ai in enumerate(atms):
        for j,aj in enumerate(atms):
            if ai>=aj:continue
            dist=mol.DM[ai-1,aj-1]
            if dist>1.7*1.889:continue
            BM[i,j]=1.0
            BM[j,i]=1.0
            vi=atomBases[ai][:,-1]
            vj=atomBases[aj][:,-1]
            if vector_angle(vi,vj)>0.5:
                print(f'原子{ai}和{aj}的法向量夹角大于90度',vi,vj)
                BM[i,j]=-1.0
                BM[j,i]=-1.0
            
    es,CM=np.linalg.eigh(BM) # 矩阵对角化
    idxs=np.argsort(-es) # 占据轨道
    es=es[idxs].copy()
    CM=CM[:,idxs].copy() # 每一列对应一个特征向量
    return BM,es,CM,idxs

# def eleMat_(mol:Mol)->np.ndarray:
#     """计算与分子轨道系数矩阵对应的电子分布矩阵"""
#     # 使用法向量可以计算每个分子的pi电子分布
#     from pywfn.maths import flib
#     obts=mol.O_obts

#     nobt=len(obts)
#     CM=mol.CM.copy()
#     nmat,nobt=CM.shape
#     NM=flib.eleMat(nmat,nobt,CM,mol.SM)*mol.oE
#     return NM

def eleMat(mol:Mol)->np.ndarray:
    """计算与分子轨道系数矩阵对应的电子分布矩阵"""
    # 使用法向量可以计算每个分子的pi电子分布
    from pywfn.maths import rlib
    obts=mol.O_obts

    nobt=len(obts)
    CM=mol.CM.copy()
    NM=rlib.ele_mat_rs(CM,mol.SM) # type: ignore
    return np.array(NM)*mol.oE

def engMat(mol:Mol,NM:np.ndarray)->np.ndarray:
    """
    电子能量分布矩阵
    """
    # obts=mol.O_obts
    obtEngs=np.array(mol.obtEngs).reshape(1,-1)
    obtOccs=np.array(mol.obtOccs).reshape(1,-1)
    EM=NM*obtEngs*obtOccs
    return EM

def piEleMat(mol:Mol)->np.ndarray:
    """计算与分子轨道系数矩阵对应的pi电子分布矩阵"""
    from pywfn.maths import rlib
    from pywfn.atomprop import direction
    dirCaler=direction.Calculator(mol)
    dirs=[]
    atms=[]
    for i in range(mol.atoms.natm):
        normal=dirCaler.normal(i+1)
        if normal is None:continue
        dirs.append(normal)
        atms.append(i+1)
    dirs=np.array(dirs)
    obts=mol.O_obts+mol.V_obts
    CMp=projCM(mol,obts,atms,dirs,False,False)
    nmat,nobt=CMp.shape
    NM=rlib.ele_mat_rs(CMp,mol.SM) # type: ignore
    return np.array(NM)*mol.oE