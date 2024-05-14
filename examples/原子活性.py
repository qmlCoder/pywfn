import sys
sys.path.append("D:\code\pywfn")

from pywfn.base import Mol
from pywfn.reader import LogReader
from pywfn.atomprop import activity

# root=rf"D:\BaiduSyncdisk\Articles\HFV\gfile\M4\M4_wfn"

# pathn=rf"{root}-.log"
# path0=rf"{root}0.log"
# pathp=rf"{root}+.log"

# mol0=Mol(reader=LogReader(path0))
# moln=Mol(reader=LogReader(pathn))
# molp=Mol(reader=LogReader(pathp))

mol=Mol(LogReader("D:\gfile\gs\int1-spin.log"))

caler=activity.Calculator()

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
# result=caler.dirFukui([2,1,3,4])
# print(result)
# for a,x,y,z,e,n in result:
#     print(f'[{a:.0f},{x:>6.2f},{y:>6.2f},{z:>6.2f},{e:>6.2f}],')
# # 计算化合价
print('自由价')
caler.mols=[mol]
result=caler.freeValence([15,18])
for a,x,y,z,v in result:
    print(f'[{a:.0f},{x:>6.2f},{y:>6.2f},{z:>6.2f},{v:>6.2f}],')

# # 计算化合价
print('自由价')
caler.mols=[mol]
result=caler.valence()
print(result[14])
print(result[17])

# print('原子能差')
# caler.mols=[moln,mol0,molp]
# engs=caler.delEng()[[36,38],:]
# for en,e0,ep in engs:
#     rn=en-e0
#     rp=e0-ep
#     print(f'{en:>6.3f},{e0:>6.3f},{ep:>6.3f},{rn:>6.3f},{rp:>6.3f}')