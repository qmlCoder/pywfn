import sys;sys.path.append("d:/code/pywfn")

import numpy as np
import pyvista as pv

from pywfn.maths import rlib
size=3
grids=[]
for x in np.linspace(-1,1,size):
    for y in np.linspace(-1,1,size):
        for z in np.linspace(-1,1,size):
            grids.append([x,y,z])


grids=np.array(grids)

vals=np.linalg.norm(grids,axis=1)
shape=[size,size,size]
verts,faces=rlib.march_cube(shape,grids.tolist(),vals.tolist(),0.5) # type: ignore
print(verts)
print(faces)
for i,xyz in enumerate(verts):
    print(i,xyz)
# pl=pv.Plotter()
# pl.add_points(grids)
# pl.add_points(np.array(verts),color="red")
# pl.show()