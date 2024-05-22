import numpy as np

from pywfn.data import sphGrid
from scipy.interpolate import griddata
def interp(gridData):
    points=gridData[:,:3]
    values=gridData[:,3]
    print(points.shape)
    new_values=griddata(points,values,new_points,method='linear')
    return new_values

new_points=np.random.rand(4,3)

gridData=sphGrid.gridData.copy()
new_values=interp(gridData)
print(new_values)

# gridData=np.load('grid_74,80.npy')
# new_values1=interp(gridData)
# print(new_values1)

# gridData=np.load('grid_74,100.npy')
# new_values2=interp(gridData)
# print(new_values2)

# print(new_values2/new_values)