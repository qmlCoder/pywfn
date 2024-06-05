import sys
sys.path.append(rf"D:\code\pywfn")

from pywfn.base import Mol
from pywfn.reader import LogReader
from pywfn.atomprop import activity
from pywfn.bondprop import bondOrder
from pywfn.atomprop import energy
import numpy as np
path0=rf"D:\gfile\ckh\M5.log"
pathn=rf"D:\gfile\ckh\M5+1.log"
pathp=rf"D:\gfile\ckh\M5-1.log"
path=rf"D:\BaiduSyncdisk\Articles\HFV\gfile\CH4\CH4_STO3.out"
path0=rf"D:\BaiduSyncdisk\Articles\HFV\gfile\cnhn\C6H6_1.log"

mol0=Mol(LogReader(path0))
# moln=Mol(LogReader(pathn))
# molp=Mol(LogReader(pathp))

# caler=activity.Calculator()

# caler.mols=[moln,mol0,molp]
# result=caler.fukui(chrg='mulliken')
# for a,atom in enumerate(mol0.atoms):
#     print(atom,result[a])

# caler.mols=[mol0]
# print(9,caler.freeValence(9))
# print(38,caler.freeValence(38))

# result=caler.dirFukui([9,38])
# print(result)
# from pywfn.atomprop import direction

# dirCaler=direction.Calculator(mol0)

# atms=[]
# dirs=[]
# for atom in mol0.atoms:
#     normal=dirCaler.normal(atom.idx) # 原子的法向量
#     if normal is None:continue
#     atms.append(atom.idx)
#     dirs.append(normal)
# CMp=mol0.projCM(mol0.O_obts,atms,dirs,False,False)

# from pywfn.molprop import energy

# eleMat=energy.Calculator(mol0).eleMat(CMp)

# print(np.sum(eleMat,axis=0))

result=energy.Calculator(mol0).atmPiEngs()
for a,atom in enumerate(mol0.atoms):
    print(atom,result[a])