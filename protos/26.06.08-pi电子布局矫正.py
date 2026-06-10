"""
计算pi电子布居的时候，需要排除不成键的分子轨道中原子的贡献
循环每个分子轨道，求轨道键级，如果原子轨道键级之和为0，则该分子轨道中该原子的贡献为0
"""

import numpy as np
from pywfn.atomprop import charge
from pywfn.base import Mole
from pywfn.bondprop import order
from pywfn.orbtprop import obtmat

mole = Mole.from_file(r"E:\code\articles\MOCV\moles\phx\phcl.fch")
mole = Mole.from_file(r"E:\code\articles\MOCV\moles\sicl2o.log")
# mole = Mole.from_file(r"C:\Users\11032\Desktop\gfile\pywfn\C6H6.out")
cmat_caler = obtmat.Calculator(mole)
stms, atos, cmat = cmat_caler.pi_mocv(0)
print(atos)
# dirs, cmat = cmat_caler.pi_pocv()
obts = mole.coefs.obt_idxs(True)  # 占据分子轨道的索引
smat = mole.basis.get_smat()  # 重叠矩阵
atm_uplo = mole.basis.atm_uplo()
obt_occs = mole.coefs.obt_occs
# cmat = mole.coefs.get_cmat()
natm = mole.atoms.len()
atm_eles = np.zeros(natm)
for obt in obts:
    cmat_part = cmat[:, obt].reshape(-1, 1)
    dmat = (cmat_part @ cmat_part.T) * obt_occs[obt]
    omat = order.mayer(dmat, smat)
    omat = order.ato2atm(omat, atm_uplo)
    weits = np.zeros(natm)
    for atm in range(natm):
        nebs = mole.atoms.nebs(atm)
        weit = omat[atm, nebs].max()
        weits[atm] = 1.0 if weit > 1e-4 else 0.0
    eles = charge.mulliken(dmat, smat, atm_uplo)
    eles = np.array(eles) * weits
    print(f"{obt + 1:>3} {eles.sum():>.2f}", weits)
    atm_eles += eles

print(atm_eles)
print(atm_eles.sum())
