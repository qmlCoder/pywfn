import sys
sys.path.append("D:\code\pywfn")
from pywfn.base import Mol
from pywfn.reader import LogReader
from pywfn.atomprop import atomCharge,atomSpin
import numpy as np
pathn=r"D:\BaiduSyncdisk\Articles\HFV\gfile\NHC\BHC_t_wfn-.log"
path0=r"D:\BaiduSyncdisk\Articles\HFV\gfile\NHC\BHC_t_wfn0.log"
pathp=r"D:\BaiduSyncdisk\Articles\HFV\gfile\NHC\BHC_t_wfn+.log"

moln=Mol(reader=LogReader(pathn))
mol0=Mol(reader=LogReader(path0))
molp=Mol(reader=LogReader(pathp))


calern=atomCharge.Calculator(moln)
caler0=atomCharge.Calculator(mol0)
calerp=atomCharge.Calculator(molp)

result=calern.mulliken()
print(result)

calern=atomSpin.Calculator(moln)
caler0=atomSpin.Calculator(mol0)
calerp=atomSpin.Calculator(molp)

result=calern.calculate()
