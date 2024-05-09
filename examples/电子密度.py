import sys

sys.path.append("D:\code\pywfn")

from pywfn.data import sphGrid
from pywfn.base import Mol
from pywfn.reader import LogReader
from pywfn.atomprop import atomCharge

import matplotlib.pyplot as plt
import pyvista as pv
import numpy as np
import time

weight = sphGrid.gridData[:, -1]
coords = sphGrid.gridData[:, :3]

path = "D:\BaiduSyncdisk\Articles\HFV\gfile\CH4.log"
# path="D:\BaiduSyncdisk\Articles\HFV\gfile\cnhn\C6H6_1.log"
# path="D:\BaiduSyncdisk\gfile\C=C\CH2=CH2.out"
path="D:\BaiduSyncdisk\gfile\elements\H2.out"
mol = Mol(reader=LogReader(path))
mol.bohr = True
atms = mol.obtAtms
shls = mol.obtShls
angs = mol.obtAngs
lmns = mol.basis.numAng(angs)

nmat = mol.CM.shape[0]
npos=len(weight)
nobt = len(mol.O_obts)
obts = mol.O_obts

wfns=np.zeros(shape=(nmat,npos))

def basDens(u):
    atm = atms[u]
    pos=mol.atom(atm).coord+coords
    wfn_u=get_wfn(u,pos)
    dens=np.zeros_like(wfn_u)
    for v in range(nmat):
        wfn_v=get_wfn(v,pos)
        dens+=wfn_u*wfn_v*weight*mol.PM[u,v]
    return np.sum(dens)


def get_wfn(i,pos):
    lmn = lmns[i]
    atm = atms[i]
    shl = shls[i]
    ang = sum(lmn)
    atmic = mol.atom(atm).atomic
    basis = mol.basis.get(atmic, shl, ang)
    exps = [b.exp for b in basis]
    coes = [b.coe for b in basis]
    pos=mol.atom(atm).coord+coords
    R2 = np.sum(pos**2, axis=1)
    wfn = mol.gto.cgf(exps, coes, lmn, R2, pos)  # 空间坐标-以原子为中心的坐标
    return wfn


def get_smi(u, v, pos):  # 重叠矩阵矩阵元
    wfn_u = get_wfn(u,pos)
    wfn_v = get_wfn(v,pos)
    smi = np.sum(wfn_u * wfn_v * weight)
    # print(u,v,smi)
    return smi


# # 重叠矩阵应该对才行
fig,axs=plt.subplots(1,2)
print(f'{nmat=}')
SM=np.zeros(shape=(nmat,nmat))
for u in range(nmat):
    for v in range(nmat):
        smi=get_smi(u,v,coords)
        SM[u,v]=smi
axs[0].imshow(SM)
axs[1].imshow(mol.SM)
plt.show()




pass
