from pywfn.base import Mol
from pywfn.maths import vector_angle

import numpy as np

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