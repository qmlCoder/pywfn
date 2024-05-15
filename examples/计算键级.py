import sys
sys.path.append("D:\code\pywfn")

from pywfn.base import Mol
from pywfn.reader import LogReader
from pywfn.atomprop import activity
from pywfn.bondprop import bondOrder


path0=fr"D:\BaiduSyncdisk\Articles\HFV\gfile\NHC\BHC_t_wfn0.log"
path0=fr"D:\BaiduSyncdisk\Articles\HFV\gfile\cnhn\C6H6_1.log"
# path0="D:\BaiduSyncdisk\gfile\scans\lianxi\dingerxi.log"
mol0=Mol(reader=LogReader(path0))

caler=bondOrder.Calculator(mol0)

result=caler.mayer()
print(result)

# result=caler.dirMayer([
#     [1,2],
#     [2,3]
# ])
# print(result)

result=caler.piOrder()
print(result)

result=caler.hmo()
idxs=(result[:,0]==1) & (result[:,1]==2)
print(result[idxs,:])
pass