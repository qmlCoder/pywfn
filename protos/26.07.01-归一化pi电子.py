from pywfn.base import Mole
from pywfn.orbtprop import obtmat

path = r"C:\Users\11032\Desktop\gfile\pywfn\C6H6.fch"
path = r"E:\code\articles\Research Diary\B12H12 Study\moles\B12H12_3n.log"
mole = Mole.from_file(path)
smat = mole.basis.get_smat()
cmat = mole.coefs.get_cmat()
caler = obtmat.Calculator(mole)
stms, atos, cmat_p = caler.pi_mocv()
idxs = mole.coefs.obt_idxs(True)  # 获取占据轨道的索引
# print(idxs, len(idxs))

cmat = cmat[:, idxs]
cmat_p = cmat_p[:, idxs]
cmat_s = cmat - cmat_p

qmat_pp = (cmat_p @ cmat_p.T) @ smat
qmat_ss = (cmat_s @ cmat_s.T) @ smat
qmat = (cmat @ cmat.T) @ smat
k = 1.0 if mole.is_open() else 2.0
eles_pp = qmat_pp.diagonal() * k
eles_ss = qmat_ss.diagonal() * k
eles = qmat.diagonal() * k
eles_p = eles_pp / (eles_pp + eles_ss) * eles
eles_s = eles_ss / (eles_pp + eles_ss) * eles

uplo = mole.basis.atm_uplo()
total = 0
for u, l in uplo:
    ele_p = eles_p[u:l].sum()
    ele_s = eles_s[u:l].sum()
    ele = eles[u:l].sum()
    total += ele_p
    print(f"{ele_p:>10.4f}, {ele_s:>10.4f}, {ele:>10.4f}")
print(total)
print("-" * 30)
total = 0
for u, l in uplo:
    ele_pp = eles_pp[u:l].sum()
    ele_ss = eles_ss[u:l].sum()
    ele = eles[u:l].sum()
    ele_p = ele_pp / (ele_pp + ele_ss) * ele
    ele_s = ele_ss / (ele_pp + ele_ss) * ele
    total += ele_p
    print(f"{ele_p:>10.4f}, {ele_s:>10.4f}, {ele:>10.4f}")
print(total)

# print(eles_pp.sum())
# print(eles_ss.sum())

# print(eles_p.sum(), eles_s.sum())
