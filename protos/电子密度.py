import sys;sys.path.append("d:/code/pywfn")
from pywfn.base import Mol
from pywfn.reader import LogReader,FchReader
from pywfn.spaceprop import density,RectGrid
import numpy as np
import matplotlib.pyplot as plt

path=rf"D:\gfile\C6H6.fch"

reader=FchReader(path)
mol=Mol(reader)

caler=density.Calculator(mol)

x,y,z=grids=np.random.rand(3)
print(f'{x},{y},{z}')
dens0,dens1,dens2=caler.molDens(grids.reshape(1,3),1)
print(dens0)
print(dens1)
print(dens2)
dens0,dens1,dens2=caler.molDens_(grids.reshape(1,3),1)
print(dens0)
print(dens1)
print(dens2)

# cent=np.array([0.,0.,0.5])
# norm=np.array([0.,0.,1.])
# vx=np.array([1,0.,0.])
# shape,grids=RectGrid().set_v1(cent,norm,vx,10).get()

# dens0,dens1,dens2=caler.molDens(grids,2)
# print(dens1.shape)
# # fig,axs=plt.subplots(1,2)
# # axs[0].imshow(dens0.reshape(shape))
# # axs[1].imshow(np.linalg.norm(dens1,axis=1).reshape(shape))
# # plt.show()

# rdg=caler.RDG(grids)
# plt.matshow(rdg.reshape(shape))
# plt.show()