import sys
sys.path.append("D:\code\pywfn")
from pywfn.base import Mol
from pywfn.reader import LogReader
from pywfn.atomprop import lowdinCharge
import numpy as np
path="D:\BaiduSyncdisk\Articles\HFV\gfile\CH4.log"
mol=Mol(reader=LogReader(path))

caler=lowdinCharge.Calculator(mol)

charges=caler.calculate()
print(charges)

# SM=np.random.randn(3,3)
# v,Q=np.linalg.eig(SM)
# V=np.diag(v)
# V_=np.linalg.inv(Q)
# SM_=Q@(V@V_)
# print(SM)
# print(SM_)