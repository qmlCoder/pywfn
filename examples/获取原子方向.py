import sys
sys.path.append("D:\code\pywfn")

import pyvista as pv

from pywfn.base import Mol
from pywfn.atomprop import direction
from pywfn.reader import GjfReader

import numpy as np
import matplotlib.pyplot as plt

def dir2arrow(atm:int,dirs:np.ndarray):
    cent=mol.atom(atm).coord
    
    for dir in dirs:
        # print(f'add_arrow {atm} {x:.2f} {y:.2f} {z:.2f} {1:.2f} 0.2 Arrow.F')
        pl.add_arrows(cent,dir)
    pl.show()

pl = pv.Plotter()
pl.add_axes_at_origin()

path="D:\BaiduSyncdisk\Articles\HFV\gfile\M4\M4_wfn.gjf"
# path="D:\BaiduSyncdisk\Articles\HFV\gfile\ok\R-C6.gjf"

mol=Mol(reader=GjfReader(path))

cal=direction.Calculator(mol)

dir2arrow(26,cal.reaAro(26))
