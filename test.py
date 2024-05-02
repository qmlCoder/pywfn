from pathlib import Path
from pywfn.tools import engCorr
import numpy as np
import matplotlib.pyplot as plt


pos1=np.random.rand(3)
pos2=np.random.rand(3)
pos3=np.random.rand(3)

# pos1=pos1/np.linalg.norm(pos1)
# pos2=pos2/np.linalg.norm(pos2)
# pos3=pos3/np.linalg.norm(pos3)

x1,y1,z1=pos1
x2,y2,z2=pos2
x3,y3,z3=pos3

A = (y2 - y1)*(z3 - z1) - (z2 - z1)*(y3 - y1)
B = (x3 - x1)*(z2 - z1) - (x2 - x1)*(z3 - z1)
C = (x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1)


print(A,B,C)

A,B,C=np.cross(pos1-pos2,pos3-pos2)

print(A,B,C)