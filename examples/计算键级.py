import sys
sys.path.append("D:\code\pywfn")

from pywfn.base import Mol
from pywfn.reader import LogReader
from pywfn.atomprop import activity
from pywfn.bondprop import bondOrder


path=fr"D:\BaiduSyncdisk\Articles\HFV\gfile\NHC\BHC_t_wfn0.log"
path=fr"D:\BaiduSyncdisk\Articles\HFV\gfile\cnhn\C6H6_1.log"
path=fr"D:\BaiduSyncdisk\gfile\scans\lianxi\lianxiScan\f01.log"
path=fr"D:\BaiduSyncdisk\gfile\CnHn\C4H4.log"
path="D:\BaiduSyncdisk\gfile\scans\lianxi\dingerxi.log"
# path=rf"D:\BaiduSyncdisk\Articles\HFV\gfile\lianben\b3lyp\logs\f13.log"
path=rf"D:\BaiduSyncdisk\Articles\HFV\gfile\NHC\NHC-wfn0.log"
mol=Mol(reader=LogReader(path))

caler=bondOrder.Calculator(mol)

result=caler.boundMayer(7)
print(result)

# result=caler.mayer()
# print(result)

# result=caler.dirMayer([
#     [1,2],
#     [2,3]
# ])
# print(result)

# result=caler.piOrder()
# print(result)

# result=caler.hmo()
# print(result)
# pass