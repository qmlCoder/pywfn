import sys;sys.path.append("d/code/pywfn")

from pywfn.base import Mole
from pywfn.reader import LogReader,FchReader,GjfReader,MdeReader,MolReader,SdfReader,XyzReader

mole=Mole(LogReader(rf"d:\gfile\各种类型文件\C6H6.out"))
print(mole.geome)
print(mole.basis)
print(mole.coefs)

mole=Mole(FchReader(rf"d:\gfile\各种类型文件\C6H6.fch"))
print(mole.geome)
print(mole.basis)
print(mole.coefs)

mole=Mole(GjfReader(rf"d:\gfile\各种类型文件\C6H6.gjf"))
print(mole.geome)

mole=Mole(MdeReader(rf"d:\gfile\各种类型文件\C6H6.molden"))
print(mole.geome)
print(mole.basis)
print(mole.coefs)

mole=Mole(MolReader(rf"d:\gfile\各种类型文件\C6H6.mol"))
print("MolReader")
print(mole.geome)

mole=Mole(SdfReader(rf"d:\gfile\各种类型文件\C6H6.sdf"))
print("SdfReader")
print(mole.geome)

mole=Mole(XyzReader(rf"d:\gfile\各种类型文件\C6H6.xyz"))
print("XyzReader")
print(mole.geome)