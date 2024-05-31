import sys
sys.path.append("D:\code\pywfn")

from pywfn.base import Mol
from pywfn.reader import GjfReader

import numpy as np

path="D:\BaiduSyncdisk\Articles\HFV\gfile\芳香性\lianbenR.gjf"
mol=Mol(GjfReader(path))

rings=[
    [10,16],
    [3,6]
]
print(mol.coords/1.889)
for ring in rings:
    idxs=[e-1 for e in ring]
    coords=mol.coords[idxs,:]
    # print(coords)
    x,y,z=np.mean(coords,axis=0)/1.889
    print(f' {"Bq":<15}{x:>14.8f}{y:>14.8f}{z:>14.8f}')