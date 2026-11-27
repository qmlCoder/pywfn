# import matplotlib.pyplot as plt
import numpy as np
from pywfn.base import Mole
from pywfn.orbtprop import obtmat

path = r"C:\Users\11032\Desktop\gfile\pywfn\C6H6.fch"
path = r"E:\code\articles\Research Diary\B12H12 Study\moles\B12H12_2n.log"
mole = Mole.from_file(path)
smat = mole.basis.get_smat()
cmat_r = mole.coefs.get_cmat()
caler = obtmat.Calculator(mole)
stms, atos, cmat_p = caler.pi_mocv()

pmat_pr = cmat_p.T @ smat @ cmat_r
pmat_pp = cmat_p.T @ smat @ cmat_p
pvec_pr = np.diagonal(pmat_pr)
pvec_pp = np.diagonal(pmat_pp)
kvec = np.divide(pvec_pr, pvec_pp, out=np.zeros_like(pvec_pr), where=(pvec_pp != 0))
kmat = np.diag(kvec)
cmat_pp = cmat_p @ kmat  # 对原始的pi分子轨道进行缩放

# print(np.diagonal(cmat_p.T @ smat @ cmat_p))
# print(np.diagonal(cmat_pp.T @ smat @ cmat_pp))


def get_eles(cmat):
    idxs = mole.coefs.obt_idxs(True)
    cmat = cmat[:, idxs]
    dmat = cmat @ cmat.T
    print(dmat.shape)
    qmat = dmat @ smat
    eles = np.diagonal(qmat)
    return eles


eles = get_eles(cmat_p)
print(eles.sum())
eles = get_eles(cmat_pp)
print(eles.sum())

idxs = mole.coefs.obt_idxs(True)

cmat_s = cmat_r - cmat_p
print(np.diagonal((cmat_p[:, idxs]).T @ smat @ cmat_p[:, idxs]).sum())
print(np.diagonal((cmat_s[:, idxs]).T @ smat @ cmat_s[:, idxs]).sum())

cmat_ss = cmat_r - cmat_pp

print(np.diagonal((cmat_pp[:, idxs]).T @ smat @ cmat_pp[:, idxs]).sum())
print(np.diagonal((cmat_ss[:, idxs]).T @ smat @ cmat_ss[:, idxs]).sum())
overlap_mat = cmat_pp.T @ smat @ cmat_ss
overlap_diag = np.diag(overlap_mat)

print("=== 同轨道π-σ正交性验证 ===")
print(f"最大重叠绝对值：{np.max(np.abs(overlap_diag)):.2e}")
print(f"前10条轨道的重叠值：{overlap_diag[:10]}")
