import sys

sys.path.append("D:\code\pywfn")

from pywfn.data import sphGrid
from pywfn.base import Mol
from pywfn.reader import LogReader
from pywfn.atomprop import atomCharge




path = "D:\BaiduSyncdisk\Articles\HFV\gfile\CH4\CH4_STO3.out"
path = "D:\BaiduSyncdisk\Articles\HFV\gfile\CH4.log"
path = r"D:\BaiduSyncdisk\Articles\HFV\gfile\NHC\BHC_t_wfn0.out"
mol = Mol(reader=LogReader(path))

caler=atomCharge.Calculator(mol)
result=caler.hirshfeld()
print(result)
print(sum(result))