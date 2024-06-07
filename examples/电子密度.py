import sys
sys.path.append(rf"D:\code\pywfn")

from pywfn.data import sphGrid,radDens
from pywfn.base import Mol
from pywfn.reader import LogReader
from pywfn.atomprop import charge
from pywfn.spaceProp import density,dftgrid,wfnfunc

import matplotlib.pyplot as plt
import numpy as np
import time

path = rf"D:\BaiduSyncdisk\Articles\HFV\gfile\CH4\CH4_STO3.out"
# path = "D:\BaiduSyncdisk\Articles\HFV\gfile\CH4\CH4.log"
path=rf"D:\BaiduSyncdisk\gfile\elements\C.out"
# path = "D:\BaiduSyncdisk\gfile\elements\O2.out"
# path = "D:\BaiduSyncdisk\gfile\elements\S2.out"
# path = "D:\BaiduSyncdisk\gfile\elements\S.out"
# path = "D:\BaiduSyncdisk\gfile\elements\CO.out"
# path="D:\BaiduSyncdisk\gfile\elements\H2.out"
# path="D:\BaiduSyncdisk\gfile\elements\He2.out"

mol = Mol(reader=LogReader(path))


wfnCaler=wfnfunc.Calculator(mol)
denCaler=density.Calculator(mol)
gridCaler=dftgrid.Calculator(mol)
# for atom in mol.atoms:
#     atm=atom.idx
#     cords,weits=dftgrid.Calculator(mol).dftGrid(atm)
#     wfns=wfnCaler.atoWfn(atm,cords)
#     print(np.sum(wfns**2*weits))

def count(array:np.ndarray,nums:list[float]):
    counts=[]
    for i in range(len(nums)-1):
        l=nums[i]
        u=nums[i+1]
        where=np.argwhere((array>l)&(array<=u))
        count=where.shape[0]
        counts.append(count)
    return counts


molGrid,molWeit=denCaler.molGrid
# print(molGrid)

# dens=denCaler.molDens_obt(molGrid)
# print(np.sum(dens*molWeit)) # 分子电子密度

# dens=denCaler.molDens_atm(molGrid)
# print(np.sum(dens*molWeit)) # 分子电子密度

# dens=denCaler.molDens_lib(molGrid.copy())
# print(np.sum(dens*molWeit)) # 分子电子密度

# grid,weit=gridCaler.dftGrid(1)

# dens=denCaler.molDens_lib(grid.copy())
# print(np.sum(dens*weit)) # 分子电子密度

grid,weit=sphGrid.grids,sphGrid.weits
rads=np.linalg.norm(grid,axis=1)
rads=np.sort(rads)
plt.plot(np.arange(len(rads)),rads)
plt.show()
print(len(molWeit))
print(count(molWeit,list(np.linspace(0,1e-4,10))))
dens=radDens.get_radDens(6,rads)


print(np.sum(dens*weit))
# cords,weits=sphGrid.cords,sphGrid.weits
# plt.scatter(np.sum(cords**2,axis=1),weits)
# plt.show()
# wfns=caler.atoWfn(1,cords)
# print(np.sum(wfns**2*weits))
# caler.molPos=cords
# wfn=caler.a2mWfn(1,cords)
# print(np.sum(wfns**2*weits))
# plt.show()
# plt.show()
pass