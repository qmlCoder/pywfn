import sys
sys.path.append("D:\code\pywfn")
from pywfn.base import Mol
from pywfn.reader import LogReader
from pywfn.atomprop import charge, spin
import numpy as np
pathn=r"D:\BaiduSyncdisk\Articles\HFV\gfile\NHC\BHC_t_wfn-.log"
path0=r"D:\BaiduSyncdisk\Articles\HFV\gfile\NHC\BHC_t_wfn0.log"
pathp=r"D:\BaiduSyncdisk\Articles\HFV\gfile\NHC\BHC_t_wfn+.log"

moln=Mol(reader=LogReader(pathn))
mol0=Mol(reader=LogReader(path0))
molp=Mol(reader=LogReader(pathp))


calern=charge.Calculator(moln)
caler0=charge.Calculator(mol0)
calerp=charge.Calculator(molp)

result=calern.mulliken()
print(result)

calern=spin.Calculator(moln)
caler0=spin.Calculator(mol0)
calerp=spin.Calculator(molp)

result=calern.calculate()
