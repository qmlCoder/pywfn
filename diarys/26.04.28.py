import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
from pywfn.base import Mole

mole = Mole.from_file(rf"c:\Users\11032\Desktop\gfile\pywfn\C6H6.out")

r0s = np.load("./r0s.npy")
rfs = np.load("./rfs.npy")
vals = np.load("./vals.npy")
cpts = np.load("./cpts.npy")

# plt.hist(vals,bins=100)
# plt.show()


# cptl=[-3,-1,1,3]
# colors=["red","green","blue","yellow"]

pl = pv.Plotter()
xyzs = mole.atoms.xyzs()
xyzs = np.array(xyzs).reshape(-1, 3)
pl.add_points(xyzs, color="black", point_size=10)
idxs = np.where(vals > 0.02)[0]
# pl.add_points(rfs[idxs], scalars=vals[idxs], cmap="coolwarm")
# pl.show()

finds = []
for idx in idxs:
    add = True
    for ref in finds:
        dist = np.linalg.norm(rfs[idx] - rfs[ref])
        if dist < 1e-6:
            add = False
            break
    if add:
        finds.append(idx.item())
print(finds)
