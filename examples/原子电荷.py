import sys

sys.path.append("D:\code\pywfn")

# from pywfn.data import sphGrid
from pywfn.base import Mol
from pywfn.reader import LogReader
from pywfn.atomprop import charge

# path = "D:\BaiduSyncdisk\Articles\HFV\gfile\CH4\CH4_STO3.out"
# path = "D:\BaiduSyncdisk\Articles\HFV\gfile\CH4.log"

root=rf"D:\BaiduSyncdisk\Articles\HFV\gfile\M4\M4_wfn"
root=rf"D:\BaiduSyncdisk\Articles\HFV\gfile\NHC\BHC_t_wfn"
root=rf"D:\BaiduSyncdisk\Articles\HFV\gfile\反极性卡宾\N2R_t_wfn"
root=rf"D:\BaiduSyncdisk\Articles\HFV\gfile\反极性卡宾\since4_s_wfn"

pathn=rf"{root}-.log"
path0=rf"{root}0.log"
pathp=rf"{root}+.log"

# path = rf"D:\BaiduSyncdisk\Articles\HFV\gfile\NHC\BHC_t_wfn0.out"

import matplotlib.pyplot as plt
results=[]
fig,axs=plt.subplots(1,2,figsize=(12,4))
cs=['b','g','r']
fs=['-','0','+']
for p,path in enumerate([pathn,path0,pathp]):
    mol = Mol(reader=LogReader(path))
    caler=charge.Calculator(mol)
    result=caler.dirCharge('mulliken',[13])[:,-1]
    results.append(result)
    print(result)
    print('-'*70)
    axs[0].plot(result,label=f'{fs[p]}',c=cs[p])
axs[0].legend()
axs[1].plot(results[0]-results[1],label='E',c='b')
axs[1].plot(results[1]-results[2],label='N',c='r')
axs[1].legend()
plt.show()

# pi电子数

# path=rf"D:\BaiduSyncdisk\Articles\HFV\gfile\反极性卡宾\N2R_t_wfn0.log"
# mol = Mol(reader=LogReader(path))
# caler=charge.Calculator(mol)
# result=caler.piElectron()
# print(result)