import sys
sys.path.append("D:\code\pywfn")

from pywfn.base import Mol
from pywfn.reader import LogReader
from pywfn.atomprop import energy
import matplotlib.pyplot as plt
import numpy as np




root='D:/BaiduSyncdisk/gfile/scans/lianxi/lianxiScan'

energList=[]
angleList=[]
energys=[]
atms=[5,1,3,8]
for i in range(36):
    path=f'{root}/f{i+1:0>2}.log'
    mol=Mol(LogReader(path))
    caler=energy.Calculator(mol)
    result=caler.dirEnergy()
    energList.append(result)

    angle=mol.params(atms)
    angleList.append(angle)
    energys.append(mol.energy)
energList=np.vstack(energList)
angleList=np.array(angleList)
fig,axs=plt.subplots(1,2,figsize=(10,6))
for atm in atms:
    axs[0].plot(angleList,energList[:,atm-1],'o-',label=f'{atm}')
axs[0].legend()
axs[1].plot(angleList,np.sum(energList,axis=1),'o-',label='total')
axs[1].legend()
plt.show()
pass