import sys
sys.path.append(rf'D:\code\pywfn')
from pywfn.base import Mol
from pywfn.reader import LogReader
from pywfn.bondprop import bondOrder

import numpy as np

root=rf"D:\BaiduSyncdisk\Articles\HFV\gfile\lianbenR\b3lyp\logs"
orderList=[]
angleList=[]
bonds=[]
for i in range(37):
    path=f'{root}/f{i+1:0>2}.log'
    mol=Mol(reader=LogReader(path))
    caler=bondOrder.Calculator(mol)
    orders=caler.piOrder()
    orderl=[]
    for a1,a2,val in orders:
        if mol.atom(int(a1)).symbol!='C':continue
        if mol.atom(int(a2)).symbol!='C':continue
        orderl.append(val)
        if i==0:bonds.append([a1,a2])
    orderList.append(orderl)
orderList=np.array(orderList)
angleList=np.array(angleList)

import matplotlib.pyplot as plt
for i in [4,7]:
    plt.plot(angleList,orderList[:,i])
plt.show()