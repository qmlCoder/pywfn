import sys;sys.path.append(rf'D:\code\pywfn') # 让解释器能找到pywfn包

from pywfn.base import Mole
from pywfn.reader import LogReader,FchReader
from pywfn.atomprop import charge,activity

reader=LogReader(rf"d:\gfile\pywfn\C6H6.out")
mole=Mole(reader)


print("原子电荷")
caler=charge.Calculator(mole)
print("mulliken")
print(caler.mulliken())
print("lowdin")
print(caler.lowdin())
# print("hirshfeld")
# print(caler.hirshfeld())
print("pi_pocv")
dirs,vals=caler.pi_pocv()
print(dirs)
print(vals)
print("pi_deco")
print(caler.pi_deco())
print("spin")
print(caler.spin())

print("原子活性")

mole_0=Mole(FchReader(rf"d:\gfile\pywfn\C6H6.fch"))
mole_n=Mole(FchReader(rf"d:\gfile\pywfn\C6H6_N.fch"))
mole_p=Mole(FchReader(rf"d:\gfile\pywfn\C6H6_P.fch"))

caler=activity.Calculator(mole_0)
print("fukui")
print(caler.fukui(mole_n,mole_p))
print("fukui_dir")
print(caler.fukui_dir(0,[0.0,0.0,1.0],mole_n,mole_p))
print("fukui_pi")
dirs,vals=caler.fukui_pi(mole_n,mole_p)
print(dirs)
print(vals)
