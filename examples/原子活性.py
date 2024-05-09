import sys
sys.path.append("D:\code\pywfn")

from pywfn.base import Mol
from pywfn.reader import LogReader
from pywfn.atomprop import activity


pathn=r"D:\gfile\ckh\N-1.log"
path0=r"D:\gfile\ckh\N.log"
pathp=r"D:\gfile\ckh\N+1.log"

moln=Mol(reader=LogReader(pathn))
mol0=Mol(reader=LogReader(path0))
molp=Mol(reader=LogReader(pathp))

caler=activity.Calculator()

# 计算福井函数
caler.mols=[moln,mol0,molp]
result=caler.fukui()
print(result)

# 计算parr函数
caler.mols=[moln,molp]
result=caler.parr()
print(result)

# 计算方向福井函数
print('方向fukui')
caler.mols=[moln,mol0,molp]
result=caler.dirFukui([37,39])
print(result)

# 计算化合价
print('化合价')
caler.mols=[mol0]
result=caler.valence()
print(37,result[36])
print(39,result[38])

# 计算自由价
print('自由价')
caler.mols=[mol0]
result=caler.freeValence([37,39])
print(result)

print('原子能差')
caler.mols=[moln,mol0,molp]
engs=caler.delEng()[[36,38],:]
for en,e0,ep in engs:
    rn=en-e0
    rp=e0-ep
    print(f'{en:>6.3f},{e0:>6.3f},{ep:>6.3f},{rn:>6.3f},{rp:>6.3f}')