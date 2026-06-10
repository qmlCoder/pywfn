"""
对pi键级进行矩阵矫正
求密度矩阵的时候使用分块乘法
"""

import numpy as np
from pywfn import atomprop
from pywfn.base import Mole
from pywfn.bondprop import order
from pywfn.orbtprop import obtmat

path = r"E:\code\articles\MOCV\moles\B12H12_21_wfn.log"
mole = Mole.from_file(path)
mat_caler = obtmat.Calculator(mole)
stms, atos, cmat = mat_caler.pi_mocv(0)
uplo = mole.basis.shls.uplo()
nshl = len(uplo)
obt_idxs = mole.coefs.obt_idxs(True)
# cmat = mole.coefs.get_cmat()  # 分子轨道系数矩阵

cmat = cmat[:, obt_idxs].copy()

nato = mole.basis.atos.len()

dmat = np.zeros((nato, nato))
keeps = np.array([1.0, 0.0, 0.0]).reshape(3, 1)
dir_caler = atomprop.direction.Calculator(mole)
for i in range(nshl):
    for j in range(nshl):
        atmi = mole.basis.shls.get(i).atm
        atmj = mole.basis.shls.get(j).atm
        ui, li = uplo[i]
        uj, lj = uplo[j]
        rcoefs_i = cmat[ui:li, :]
        rcoefs_j = cmat[uj:lj, :]
        if li - ui == 3 and lj - uj == 3 and atmi != atmj:
            stm = dir_caler.LCS(atmi, atmj)
            tcoefs_i = obtmat.mocv(stm.mat, rcoefs_i, keeps)
            tcoefs_j = obtmat.mocv(stm.mat, rcoefs_j, keeps)
            dmat[ui:li, uj:lj] = (tcoefs_i @ tcoefs_j.T) * 2
        else:
            dmat[ui:li, uj:lj] = (rcoefs_i @ rcoefs_j.T) * 2
            # pass
# dmat = mole.coefs.get_dmat()
smat = mole.basis.get_smat()
ato_omat = order.mayer(dmat, smat)
atm_uplo = mole.basis.atm_uplo()
atm_omat = order.ato2atm(ato_omat, atm_uplo)
natm = mole.atoms.len()
for i in range(natm):
    for j in range(natm):
        if not mole.atoms.is_bond(i, j):
            continue
        if i == j:
            continue
        val = atm_omat[i, j] ** 0.5
        print(f"atm_omat[{i}, {j}] = {val}")
