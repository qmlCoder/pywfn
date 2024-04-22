import sys
sys.path.append("D:\code\pywfn")

from pywfn.base import Mol
from pywfn.reader import LogReader

path="D:\BaiduSyncdisk\Articles\HFV\gfile\CH4.log"
reader=LogReader(path)

# occs=reader.get_obtOccs()
# print(occs)
# engs=reader.get_obtEngs()
# print(engs)
# atms=reader.get_obtAtms()
# print(atms)
# angs=reader.get_obtAngs()
# print(angs)
basis=reader.read_basisData()
for each in basis:
    print(each)

mol=Mol(reader=reader)
mol.gto.bind(1,1)