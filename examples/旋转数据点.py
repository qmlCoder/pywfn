import sys
sys.path.append("D:\code\pywfn")

from pywfn import maths
import numpy as np
import pyvista as pv

pl = pv.Plotter()
pl.add_axes_at_origin()



center=np.random.rand(3)

points=np.random.rand(1,3)+center

anyv=np.random.rand(3)+center

axis=np.cross(anyv,points.flatten()-center)
angle=0.1

for ang in np.linspace(0,np.pi*2,10,endpoint=False):
    points2=maths.points_rotate(points,center,axis,ang)
    pl.add_arrows(center.reshape(1,3),points2-center)

pl.show()