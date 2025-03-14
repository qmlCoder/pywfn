import sys;sys.path.append("d:/code/pywfn")
from pywfn.base import Mol
from pywfn.reader import LogReader,FchReader
from pywfn.spaceprop import density,RectGrid,dftgrid
import numpy as np
import matplotlib.pyplot as plt

path=rf"D:\gfile\H2.fch"
path=rf"D:\gfile\H2O_b3lyp.fch"

reader=FchReader(path)
mol=Mol(reader)

caler=dftgrid.Calculator(mol)
# caler.nrad=1
# caler.nsph=6

res1=caler.a2mGrid(1)[1]

res2=caler.a2mGrid_(1)[1]
print(res1)
print(res2)