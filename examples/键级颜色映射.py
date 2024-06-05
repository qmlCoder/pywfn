import sys
sys.path.append("D:\code\pywfn")

from pywfn.base import Mol
from pywfn.reader import LogReader
from pywfn.atomprop import activity
from pywfn.bondprop import bondOrder

import numpy as np


path=fr"D:\BaiduSyncdisk\Articles\HFV\gfile\NHC\BHC_t_wfn0.log"
path=fr"D:\BaiduSyncdisk\Articles\HFV\gfile\cnhn\C6H6_1.log"
path=fr"D:\BaiduSyncdisk\gfile\scans\lianxi\lianxiScan\f01.log"
path=fr"D:\BaiduSyncdisk\gfile\CnHn\C4H4.log"
path="D:\BaiduSyncdisk\gfile\scans\lianxi\dingerxi.log"
# path=rf"D:\BaiduSyncdisk\Articles\HFV\gfile\lianben\b3lyp\logs\f13.log"
path=rf"D:\BaiduSyncdisk\Articles\HFV\gfile\NHC\NHC-wfn0.log"
path=rf"D:\BaiduSyncdisk\gfile\scans\lianbenR\lianbenR_wfn.log"
# path=rf"D:\BaiduSyncdisk\gfile\CnHn\C6H6.log"
mol=Mol(reader=LogReader(path))

caler=bondOrder.Calculator(mol)

# result=caler.boundMayer(7)
# print(result)

# result=caler.mayer()
# print(result)

# result=caler.dirMayer([
#     [1,2],
#     [2,3]
# ])
# print(result)

def cmap(val):
    c0=np.array([0,0,1])
    c1=np.array([1,0,0])
    dc=c1-c0
    return c0+dc*val

result=caler.piOrder()
# print(result)
vals=[]
for a1,a2,val in result:
    if mol.atom(int(a1)).symbol=='H':continue
    if mol.atom(int(a2)).symbol=='H':continue
    vals.append(val)
    r,g,b=cmap(val)*255
    a1=int(a1)
    a2=int(a2)
    print(f'{a1:>2}-{a2:>2}:{val:.2f}    {r:>6.2f},{g:>6.2f},{b:>6.2f}')
print(np.std(vals))

# result=caler.hmo()
# print(result)
# pass