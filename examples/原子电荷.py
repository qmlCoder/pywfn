import sys
sys.path.append("D:\code\pywfn")
from pywfn.base import Mol
from pywfn.reader import LogReader
from pywfn.atomprop import atomCharge
import numpy as np
path="D:\BaiduSyncdisk\Articles\HFV\gfile\CH4.log"
path="D:\BaiduSyncdisk\Articles\HFV\gfile\cnhn\C6H6_1.log"
mol=Mol(reader=LogReader(path))

caler=atomCharge.Calculator(mol)

charges=caler.hirshfeld()
print(charges)

# SM=np.random.randn(3,3)
# v,Q=np.linalg.eig(SM)
# V=np.diag(v)
# V_=np.linalg.inv(Q)
# SM_=Q@(V@V_)
# print(SM)
# print(SM_)