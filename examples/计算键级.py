import sys
sys.path.append("D:\code\pywfn")

from pywfn.base import Mol
from pywfn.reader import LogReader
from pywfn.atomprop import activity
from pywfn.bondprop import bondOrder


path0=r"D:\BaiduSyncdisk\Articles\HFV\gfile\NHC\BHC_t_wfn0.log"
mol0=Mol(reader=LogReader(path0))

caler=bondOrder.Calculator(mol0)

result=caler.mayer()
print(result)

result=caler.dirMayer([
    [2,3],
    [1,6]
])
print(result)