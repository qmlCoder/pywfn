"""
直接使用变换矩阵的话，球谐的可以转换为笛卡尔的，但是笛卡尔的不能转换为球谐的，因为有冗余
"""
from pywfn.base import Mole
from pywfn.datas import btrans
import numpy as np

Pi = np.array([
    [2/3, -1/3, -1/3, 0, 0, 0],
    [-1/3, 2/3, -1/3, 0, 0, 0],
    [-1/3, -1/3, 2/3, 0, 0, 0],
    [0,0,0,1,0,0],
    [0,0,0,0,1,0],
    [0,0,0,0,0,1]
])

print(Pi)
print(btrans.DMAT.T@Pi)
print(btrans.DMAT@btrans.DMAT.T)

path=r"c:\Users\11032\Desktop\gfile\pywfn\C6H6.out"
mole=Mole.from_file(path)
idxs=mole.coefs.obt_idxs(True)
cmat_raw=mole.coefs.get_cmat("raw")[:,idxs]
smat_raw=mole.basis.get_smat("raw")
uplo_raw=mole.basis.shls.uplo("raw")
# print(uplo_raw)
cmat_sph=[]
for up,lo in uplo_raw:
    part=cmat_raw[up:lo,:]
    if (lo-up)==1:
        cmat_sph.append(part)
    if (lo-up)==3:
        cmat_sph.append(part)
    if (lo-up)==6:
        # print((btrans.DMAT.T@Pi)@part)
        # print(btrans.DMAT.T@part)
        cmat_sph.append((btrans.DMAT.T@Pi)@part)
        # cmat_sph.append(btrans.DMAT.T@part)
cmat_sph=np.concat(cmat_sph,axis=0)
smat_sph=mole.basis.get_smat("sph")
print(np.trace(cmat_raw.T@smat_raw@cmat_raw))
print(np.trace(cmat_sph.T@smat_sph@cmat_sph))
