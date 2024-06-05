import sys
sys.path.append("D:\code\pywfn")
from pywfn.base import Mol
from pywfn.reader import LogReader
from pywfn.atomprop import delProps,piProps,dirProps
from pywfn.molprop import energy
import numpy as np

# root=r"D:\BaiduSyncdisk\Articles\HFV\gfile\NHC\NHC-wfn"
root=r"D:\BaiduSyncdisk\Articles\HFV\gfile\M4\M4_wfn"
paths=[f'{root}{e}.log' for e in ('-','','+')]
mols=[Mol(reader=LogReader(path)) for path in paths]

caler=delProps.Calculator(mols)

caler.atoms=mols[0].atoms.indexs
atoms=[3,25,26,27]
caler.atoms=atoms
vect=np.array([0.30,0.48,0.83])
vect=vect/np.linalg.norm(vect)
caler.vects=[vect]*len(caler.atoms)
engs=caler.atomEngs()
# print(caler.allEngs)
chrgs=np.zeros(shape=(3,4)) #四个原子，三种电荷
for m,mol in enumerate(mols):
    
    caler2=dirProps.Calculator(mol)
    caler2.prop='charge'
    caler2.chrg='mulliken'
    caler2.vects=[vect]*len(caler.atoms)
    caler2.atoms=atoms
    chrgs[m]=caler2.calculate()
for e,eng in enumerate(engs):
    if e+1 not in caler.atoms:continue
    print(e+1,chrgs[:,e],caler.allEngs[e],eng)
    

# caler=obtEnergy.Calculator(mols[2])
# EM=caler.calculate()
# NM=caler.NM
# print(NM.sum(axis=0))
# print(EM.sum(axis=0))