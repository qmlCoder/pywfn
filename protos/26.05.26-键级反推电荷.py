"""
根据pi键级反推出电子布局
mayer键级公式复杂无法反推，只能根据wiberg键级反推

结论：键级反推原子电荷行不通，平方再开根号丢失正负信息!!
"""

import numpy as np
from pywfn import maths
from pywfn.base import Mole
from pywfn.bondprop import order

path = r"E:\code\articles\MOCV\moles\phx\phcl.out"
mole = Mole.from_file(path)
smat = mole.basis.get_smat("raw")
dmat = mole.coefs.get_dmat("raw")
cmat = mole.coefs.get_cmat("raw")

omat = order.wiberg(dmat, smat)

pxmat = maths.matrix_phalf(smat)
pmat = pxmat @ dmat @ pxmat
print(np.min(pmat), np.max(pmat))
pmat = pmat**2
print(np.allclose(omat, pmat))

nxmat = maths.matrix_nhalf(smat)
res = nxmat @ (pmat**0.5) @ nxmat
print(np.allclose(res, dmat))
