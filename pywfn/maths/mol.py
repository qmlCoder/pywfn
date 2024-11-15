from pywfn.base import Mol
from pywfn.maths import vector_angle
from pywfn.maths import flib

import numpy as np
import re
from itertools import product

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
        syms=mol.obtSyms[u:l] #该原子的轨道符号
        
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
            # if len(Cop)>6:Cop[:3]=0 # 只要价层轨道
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
    # 1.建立键连矩阵
    atms=mol.heavyAtoms # 重原子列表
    natm=len(atms) # 重原子数量
    BM=np.zeros(shape=(natm,natm)) # 键连矩阵
    DM=np.zeros_like(BM) # 键长矩阵
    for i,j in product(range(natm),range(natm)):
        a1,a2=atms[i],atms[j]
        if a1>=a2:continue
        dist=mol.DM[i,j]
        if dist>1.7*1.889:continue
        BM[i,j]=1.0
        BM[j,i]=1.0
    # 2.求解
    es,CM=np.linalg.eig(BM) # 矩阵对角化
    nele=int(len(atms)-mol.charge) #电子数量
    # idxs=np.argsort(es)[:nele//2] # 占据轨道
    idxs=np.argsort(-es) # 占据轨道
    es=es[idxs].copy()
    CM=CM[:,idxs].copy() # 每一列对应一个特征向量
    # sC=C[:,idxs].copy() # 每一列对应一个特征向量
    # se=e[idxs].copy()
    return BM,es,CM,idxs

def eleMat(mol:Mol)->np.ndarray:
    """计算与分子轨道系数矩阵对应的电子分布矩阵"""
    # 使用法向量可以计算每个分子的pi电子分布
    
    obts=mol.O_obts

    nobt=len(obts)
    CM=mol.CM[:,obts].copy()
    nmat,nobt=CM.shape
    NM=flib.eleMat(nmat,nobt,CM,mol.SM)*mol.oE
    return NM



def engMat(mol:Mol,NM:np.ndarray)->np.ndarray:
    """
    电子能量分布矩阵
    """
    obts=mol.O_obts
    obtEngs=np.array(mol.obtEngs).reshape(1,-1)[:,obts]
    EM=NM*obtEngs
    return EM

def piEleMat(mol:Mol)->np.ndarray:
    """计算与分子轨道系数矩阵对应的pi电子分布矩阵"""
    from pywfn.atomprop import direction
    dirCaler=direction.Calculator(mol)
    dirs=[]
    atms=[]
    for i in range(mol.atoms.natm):
        normal=dirCaler.normal(i+1)
        if normal is None:continue
        dirs.append(normal)
        atms.append(i+1)
    obts=mol.O_obts+mol.V_obts
    CMp=projCM(mol,obts,atms,dirs,False,False)
    nmat,nobt=CMp.shape
    NM=flib.eleMat(nmat,nobt,CMp,mol.SM)*mol.oE
    return NM