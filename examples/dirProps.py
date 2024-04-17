import sys
sys.path.append("D:\code\pywfn")
from pywfn.base import Mol
from pywfn.reader import LogReader
from pywfn.atomprop import delProps
from pywfn.molprop import obtEnergy
import numpy as np

# root=r"D:\BaiduSyncdisk\Articles\HFV\gfile\NHC\NHC-wfn"
root=r"D:\BaiduSyncdisk\Articles\HFV\gfile\M4\M4_wfn"
paths=[f'{root}{e}.log' for e in ('-','','+')]
mols=[Mol(reader=LogReader(path)) for path in paths]

caler=delProps.Calculator(mols)

caler.atoms=mols[0].atoms.indexs
caler.atoms=[3,25,26,27]
vect=np.array([0.30,0.48,0.83])
vect=vect/np.linalg.norm(vect)
# caler.vects=[vect]*len(caler.atoms)
engs=caler.atomEngs()
print(caler.allEngs)
for e,eng in enumerate(engs):
    print(e+1,eng)

# caler=obtEnergy.Calculator(mols[2])
# EM=caler.calculate()
# NM=caler.NM
# print(NM.sum(axis=0))
# print(EM.sum(axis=0))