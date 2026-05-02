import pyvista as pv
import numpy as np
from pywfn.base import Mole
import matplotlib.pyplot as plt

mole=Mole.from_file(rf"c:\Users\11032\Desktop\gfile\pywfn\C6H6.out")

r0s=np.load("./r0s.npy")
rfs=np.load("./rfs.npy")
vals=np.load("./vals.npy")
cpts=np.load("./cpts.npy")

# plt.hist(vals,bins=100)
# plt.show()


# cptl=[-3,-1,1,3]
# colors=["red","green","blue","yellow"]

pl=pv.Plotter()
xyzs=mole.atoms.xyzs()
xyzs=np.array(xyzs).reshape(-1,3)
pl.add_points(xyzs,color="black",point_size=10)
idxs=np.where(vals>0.02)[0]
pl.add_points(rfs[idxs], scalars=cpts[idxs], cmap="coolwarm")
# for i in range(len(cptl)):
#     idxs=np.where((cpts==cptl[i])&(vals>0.02))[0]
#     print(cptl[i],rfs[idxs])
    
    
#     pl.add_points(rfs[idxs],color="blue")

pl.show()
