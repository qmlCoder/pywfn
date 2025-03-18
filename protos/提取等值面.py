import sys;sys.path.append("d:/code/pywfn")

import numpy as np
import pyvista as pv


from pywfn.maths import rlib
from pywfn.base import Mole
from pywfn.reader import LogReader
from pywfn.gridprop import density,CubeGrid

size=30
grids=[]
vals=[]
for x in np.linspace(-1,1,size):
    for y in np.linspace(-1,1,size):
        for z in np.linspace(-1,1,size):
            grids.append([x,y,z])
            if z==0:
                vals.append(0)
            else:
                vals.append((x**2+y**2+z**2)*z/abs(z))


grids=np.array(grids)
vals=np.array(vals)
# vals=np.linalg.norm(grids,axis=1)
print(vals)
shape=[size,size,size]

# path=rf"D:\gfile\C6H6.out"
# mol=Mole(LogReader(path))
# caler=density.Calculator(mol)

# p0,p1=mol.spaceBorder
# print(p0,p1)
# shape,grids=CubeGrid().set_v1(p0,p1,0.3,2).get()

# vals=caler.molDens(grids,level=0)[0]

pverts,pfaces,nverts,nfaces=rlib.march_cube_rs(shape,grids,vals,0.5) # type: ignore

pl=pv.Plotter()
# pl.add_points(grids)
pl.add_points(np.array(pverts),color="red")
pl.add_points(np.array(nverts),color="blue")
pl.show()