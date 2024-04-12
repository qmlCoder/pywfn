from pywfn.base import Mol
from pywfn.reader import LogReader
from pywfn.atomprop import piProps,dirProps
from pywfn.bondprop import piDM
import numpy as np

path="D:\BaiduSyncdisk\Articles\HFV\gfile\cnhn\C6H6_1.log"
mol=Mol(reader=LogReader(path))
caler=piProps.Calculator(mol)
caler.chrg='mulliken'
caler.prop='charge'
res=caler.calculate()
print(res)

caler=dirProps.Calculator(mol)
caler.atoms=[1]
caler.vects=[np.array([0.,0.,1.])]
res=caler.projVector()
print(res)

caler=piDM.Calculator(mol)
caler.bond=[1,2]
res=caler.calculate()
print(res)
