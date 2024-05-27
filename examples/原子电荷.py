import sys

sys.path.append("D:\code\pywfn")

# from pywfn.data import sphGrid
from pywfn.base import Mol
from pywfn.reader import LogReader
from pywfn.atomprop import charge

# path = "D:\BaiduSyncdisk\Articles\HFV\gfile\CH4\CH4_STO3.out"
path = "D:\BaiduSyncdisk\Articles\HFV\gfile\CH4\CH4.log"
path="D:\BaiduSyncdisk\Articles\HFV\gfile\CH4\CH4_STO3.out"
# path="D:\BaiduSyncdisk\Articles\HFV\gfile\cnhn\C6H6_1.log"

mol=Mol(LogReader(path))

caler=charge.Calculator(mol)

result=caler.mulliken()
print(result)

result=caler.lowdin()
print(result)

# result=caler.sapce()
# print(result)

result=caler.hirshfeld()
print(result)
print(sum(result))