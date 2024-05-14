import sys

sys.path.append("D:\code\pywfn")

from pywfn.data import sphGrid
from pywfn.base import Mol
from pywfn.reader import LogReader
from pywfn.atomprop import charge

import matplotlib.pyplot as plt
import pyvista as pv
import numpy as np
import time

weight = sphGrid.gridData[:, -1]
coords = sphGrid.gridData[:, :3]

path = "D:\BaiduSyncdisk\Articles\HFV\gfile\CH4\CH4_STO3.out"
path = "D:\BaiduSyncdisk\Articles\HFV\gfile\CH4.log"
# path="D:\BaiduSyncdisk\Articles\HFV\gfile\cnhn\C6H6_1.log"
# path="D:\BaiduSyncdisk\gfile\C=C\CH2=CH2.out"
# path = "D:\BaiduSyncdisk\gfile\elements\H2.out"
mol = Mol(reader=LogReader(path))

atms = mol.obtAtms
shls = mol.obtShls
syms = mol.obtSyms
lmns = [mol.basis.sym2lmn(sym) for sym in syms]

nmat = mol.CM.shape[0]
npos = len(weight)
nobt = len(mol.O_obts)
obts = mol.O_obts

def basDens(u):
    atm = atms[u]
    pos = mol.atom(atm).coord + coords
    wfn_u = get_wfn(u, pos)
    dens = np.zeros_like(wfn_u)
    for v in range(nmat):
        wfn_v = get_wfn(v, pos)
        dens += wfn_u * wfn_v * weight * mol.PM[u, v]
    return np.sum(dens)


def get_wfn(i, pos):
    lmn = lmns[i]
    atm = atms[i]
    shl = shls[i]
    ang = sum(lmn)
    atom = mol.atom(atm)
    atmic = atom.atomic
    basis = mol.basis.get(atmic, shl, ang)
    exps = [b.exp for b in basis]
    coes = [b.coe for b in basis]
    pos_ = pos - atom.coord
    R2 = np.sum(pos_**2, axis=1)
    wfn = mol.gto.cgf(exps, coes, lmn, R2, pos_)  # 空间坐标-以原子为中心的坐标
    return wfn

errMat=np.zeros((nmat,nmat))
def get_smi(u, v, pos):  # 重叠矩阵矩阵元，与轨道组合系数无关
    wfn_u = get_wfn(u, pos)
    wfn_v = get_wfn(v, pos)
    
    smi = np.sum(wfn_u * wfn_v * weight)
    if abs(smi-mol.SM[u,v])>0.01:
        # print(u,v,atms[u],atms[v],np.all(wfn_u==wfn_v),f'{smi:.4f}',f'{mol.SM[u,v]:>.4f}')
        # print(shls[u],syms[u],shls[v],syms[v])
        errMat[u,v]=1
    # print(u,v,smi)
    return smi


# 分子电子密度
# molDens=np.zeros(npos)
# for o in obts:
#     obtDens=np.zeros(npos) # 每个分子轨道的密度
#     for u in range(nmat):
#         pos=coords #+mol.atom(atms[u]).coord # v和u应该使用相同的空间坐标啊
#         # pos=coords
#         wfn_u=get_wfn(u,pos)
#         for v in range(nmat):
#             wfn_v=get_wfn(v,pos)
#             obtDens+=wfn_u*wfn_v*mol.CM[u,o]*mol.CM[v,o]*2
#     print(o,np.sum(obtDens*weight))
#     molDens+=obtDens
# np.sum(molDens*weight)

## 原子电子密度
for atom in mol.atoms:
    wfn=mol.gto.ato(coords-atom.coord,atom.idx,mol.O_obts)
    print(np.sum(wfn*wfn*weight))

# 重叠矩阵应该对才行
# fig, axs = plt.subplots(1, 2)
# print(f"{nmat=}")
# SM = np.zeros(shape=(nmat, nmat))
# for u in range(nmat):
#     up = mol.atom(atms[u]).coord
#     for v in range(nmat):
#         vp = mol.atom(atms[v]).coord
#         pos = coords +(up+vp)/2
#         smi = get_smi(u, v, pos)
#         SM[u, v] = smi
# # axs[0].imshow(SM)
# # axs[1].imshow(mol.SM)

# dSM = mol.SM - SM
# print(np.min(dSM), np.max(dSM))
# plt.imshow(dSM,vmin=-1,vmax=1,cmap='hot')
# plt.colorbar()
# plt.show()

pass
