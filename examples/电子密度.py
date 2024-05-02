import sys
sys.path.append("D:\code\pywfn")

from pywfn.data import sphGrid
from pywfn.base import Mol
from pywfn.reader import LogReader
from pywfn.atomprop import atomCharge

import pyvista as pv
import matplotlib.pyplot as plt
import numpy as np

weight=sphGrid.gridData[:,-1]
coords=sphGrid.gridData[:,:3]


R2=np.sum(coords**2,axis=1)
path="D:\BaiduSyncdisk\Articles\HFV\gfile\CH4.log"
path="D:\BaiduSyncdisk\Articles\HFV\gfile\cnhn\C6H6_1.log"
mol=Mol(reader=LogReader(path))
bdata=mol.basis.get(1)
# print(bdata)

atms=mol.obtAtms
shls=mol.obtShls
angs=mol.obtAngs
lmns=mol.basis.numAng(angs)

nmat=mol.CM.shape[0]
nobt=len(mol.O_obts)
obts=mol.O_obts

for o in obts:
    dens=np.zeros(len(weight))
    for u in range(nmat):
        lmn_u=lmns[u]
        basis_u=mol.basis.get(atms[u],shls[u],sum(lmn_u))
        exps_u=[b.exp for b in basis_u]
        coes_u=[b.coe for b in basis_u]
        wfn_u=mol.gto.cgf(exps_u,coes_u,lmn_u,R2,coords)
        for v in range(nmat):
            lmn_v=lmns[v]
            basis_v=mol.basis.get(atms[v],shls[v],sum(lmn_v))
            exps_v=[b.exp for b in basis_v]
            coes_v=[b.coe for b in basis_v]
            wfn_v=mol.gto.cgf(exps_v,coes_v,lmn_v,R2,coords)
            dens+=mol.CM[u,o]*mol.CM[v,o]*wfn_u*wfn_v*2

    print(np.sum(dens*weight)) #每个分子轨道应该有两个电子才对



# params=mol.gto.bind(1,0)
# for C,exp,coe,l,m,n in params:

#     if abs(C)<1e-6:continue #系数很小的忽略
#     wfn=C*mol.gto.gto(exp,coe,R2,coords,l,m,n)
#     break
# dens=wfn**2
# plt.scatter(np.sqrt(R2),dens)
# print(np.sum(wfn**2))
# print(mol.SM[0,0])

# for obt in mol.O_obts: #每一个轨道所有原子电子密度的加和，理论上每个轨道中的电子应该为2
#     Qsum=0
#     for atom in mol.atoms:
#         dens=atom.get_dens([1],coords)
#         Q=np.sum(dens*weight)
#         Qsum+=Q
#     print(obt,Qsum)