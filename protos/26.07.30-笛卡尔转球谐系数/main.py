"""
笛卡尔系数到球谐系数的转换
https://chat.deepseek.com/share/sii80nt73tjgbjggzo
"""

from pywfn.datas.btrans import DMAT
from pywfn.base import Mole
import numpy as np

path=r"C:\Users\11032\Desktop\gfile\pywfn\C6H6.out"

mole=Mole.from_file(path)
idxs=mole.coefs.obt_idxs(True)
cmat_car=mole.coefs.get_cmat("raw")[:,idxs]
smat_car=mole.basis.get_smat("car")
smat_sph=mole.basis.get_smat("sph")

uplo_car=mole.basis.shls.uplo("car")
uplo_sph=mole.basis.shls.uplo("sph")

nsph=mole.basis.atos.len("sph")
ncar=mole.basis.atos.len("car")
nshl=mole.basis.shls.len()
T=np.zeros((ncar,nsph))

parts=[]
for i in range(nshl):
    u_car,l_car=uplo_car[i]
    u_sph,l_sph=uplo_sph[i]
    n_car=l_car-u_car
    n_sph=l_sph-u_sph

    part_car=cmat_car[u_car:l_car,:]
    if n_car==n_sph:
        T[u_car:l_car,u_sph:l_sph]=np.eye(n_car)
    else:
        if n_car==6:
            T[u_car:l_car,u_sph:l_sph]=DMAT


cmat_sph = np.linalg.inv(T.T @ smat_car @ T)@T.T@smat_car@cmat_car
# cmat_sph = np.linalg.solve(T.T @ smat_car @ T, T.T @ smat_car @ cmat_car)
print(np.diagonal(cmat_car.T @ smat_car @ cmat_car).sum())
print(np.diagonal(cmat_sph.T @ smat_sph @ cmat_sph).sum())
