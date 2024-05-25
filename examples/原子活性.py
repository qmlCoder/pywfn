import sys
sys.path.append("D:\code\pywfn")

from pywfn.base import Mol
from pywfn.reader import LogReader
from pywfn.atomprop import activity

root=rf"D:\BaiduSyncdisk\Articles\HFV\gfile\M4\M4_wfn"
root=rf"D:\BaiduSyncdisk\Articles\HFV\gfile\NHC\BHC_t_wfn"
root=rf"D:\BaiduSyncdisk\Articles\HFV\gfile\反极性卡宾\N2R_t_wfn"
# root=rf"D:\BaiduSyncdisk\Articles\HFV\gfile\反极性卡宾\since4_s_wfn"

pathn=rf"{root}-.log"
path0=rf"{root}0.log"
pathp=rf"{root}+.log"

mol0=Mol(reader=LogReader(path0))
moln=Mol(reader=LogReader(pathn))
molp=Mol(reader=LogReader(pathp))

caler=activity.Calculator()

caler.mols=[mol0]
result=caler.freeValence(1)
print(result)


# results=[]
# for mol in [moln,mol0,molp]:
#     caler.mols=[mol]
#     result=caler.freeValence([3,25,26,27])
#     results.append(result)
#     print(mol)



# mol=Mol(LogReader("D:\gfile\gs\int1-spin.log"))

# caler=activity.Calculator()

# # 计算福井函数
# caler.mols=[moln,mol0,molp]
# result=caler.fukui()
# print(result)

# # 计算parr函数
# caler.mols=[moln,molp]
# result=caler.parr()
# print(result)

# 计算方向福井函数
# print('方向fukui')
# caler.mols=[moln,mol0,molp]
# result=caler.dirFukui([13])
# print(result)

# me=min(result[:,4])
# mn=min(result[:,5])
# me=0
# mn=0
# print('亲电性')
# for a,x,y,z,e,n in result:
#     print(f'[{a:.0f},{x:>6.2f},{y:>6.2f},{z:>6.2f},{e-me:>6.2f}],')
    # print(f'[{a:.0f},{x:>6.2f},{y:>6.2f},{z:>6.2f},{n-mn:>6.2f}],')

# print('亲核性')
# for a,x,y,z,e,n in result:
#     # print(f'[{a:.0f},{x:>6.2f},{y:>6.2f},{z:>6.2f},{e-me:>6.2f}],')
#     print(f'[{a:.0f},{x:>6.2f},{y:>6.2f},{z:>6.2f},{n-mn:>6.2f}],')

# import matplotlib.pyplot as plt
# plt.plot(result[:,4],c='b')
# plt.plot(result[:,5],c='r')
# plt.show()
# # 计算化合价

# print('自由价')
# caler.mols=[mol]
# result=caler.freeValence([15,18])
# for a,x,y,z,v in result:
#     print(f'[{a:.0f},{x:>6.2f},{y:>6.2f},{z:>6.2f},{v:>6.2f}],')

# # # 计算化合价
# print('自由价')
# caler.mols=[mol]
# result=caler.valence()
# print(result[14])
# print(result[17])

# print('原子能差')
# caler.mols=[moln,mol0,molp]
# engs=caler.delEng()[[36,38],:]
# for en,e0,ep in engs:
#     rn=en-e0
#     rp=e0-ep
#     print(f'{en:>6.3f},{e0:>6.3f},{ep:>6.3f},{rn:>6.3f},{rp:>6.3f}')