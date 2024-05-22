import sys
sys.path.append("D:\code\pywfn")

from pywfn.data import sphGrid
from pywfn.base import Mol
from pywfn.reader import LogReader
from pywfn.atomprop import charge
from pywfn.spaceProp import density

import matplotlib.pyplot as plt
import pyvista as pv
import numpy as np
import time

path = "D:\BaiduSyncdisk\Articles\HFV\gfile\CH4\CH4_STO3.out"
path = "D:\BaiduSyncdisk\Articles\HFV\gfile\CH4.log"

mol = Mol(reader=LogReader(path))

caler=density.Calculator(mol)


for i in range(4):
    result=caler.atmDens(i+1,caler.points)
    print(np.sum(result*caler.weights))
pass
