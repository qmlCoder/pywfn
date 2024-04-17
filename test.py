from pywfn.reader import LogReader
from pywfn.base import Mol
from pywfn.atomprop import atomCharge,delProps
from pywfn.data import sphGrid
import numpy as np
import matplotlib.pyplot as plt

coord=sphGrid.gridData[:,:3]
length=np.linalg.norm(coord,axis=1)
weight=sphGrid.gridData[:,3]
plt.scatter(length,weight)
# print(np.sum(densi))
# for dis in dista:
    # print(dis)
# fig,axs=plt.subplots(1,3)
# x=coord[:,0]
# y=coord[:,1]
# z=coord[:,2]
# axs[0].scatter(x,y)
# axs[1].scatter(x,z)
# axs[2].scatter(y,z)
plt.show()