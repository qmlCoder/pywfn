import sys
sys.path.append("D:\code\pywfn\src")

from pywfn.base import Mol
from pywfn.reader import LogReader
from pywfn.molprop import aromaticity

path=fr"D:\BaiduSyncdisk\Articles\HFV\gfile\cnhn\C6H6_1.log"
# path=fr"D:\BaiduSyncdisk\gfile\scans\lianxi\lianxiScan\f01.log"
mol=Mol(reader=LogReader(path))

caler=aromaticity.Calculator(mol)

result=caler.calculate()

print(result)