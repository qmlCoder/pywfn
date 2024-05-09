import sys
sys.path.append("D:\code\pywfn")

from pywfn.base import Mol
from pywfn.reader import LogReader
from pywfn.atomprop import activity
from pywfn.bondprop import bondDirect

path0=r"D:\BaiduSyncdisk\Articles\HFV\gfile\NHC\BHC_t_wfn0.log"
mol0=Mol(reader=LogReader(path0))

caler=bondDirect.Calculator(mol0)

result=caler.verticals(2,3)
print(result)