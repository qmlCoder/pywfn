import sys
sys.path.append("D:\code\pywfn")

from pywfn.data import sphGrid
from pywfn.base import Mol
from pywfn.reader import LogReader
from pywfn.atomprop import charge
from pywfn.spaceProp import density

import matplotlib.pyplot as plt
import pyvista as pv
import numpy as np
import time

path = "D:\BaiduSyncdisk\Articles\HFV\gfile\CH4\CH4_STO3.out"
# path = "D:\BaiduSyncdisk\Articles\HFV\gfile\CH4\CH4.log"
# path="D:\BaiduSyncdisk\gfile\elements\C.out"
# path = "D:\BaiduSyncdisk\gfile\elements\O2.out"
# path = "D:\BaiduSyncdisk\gfile\elements\S2.out"
# path = "D:\BaiduSyncdisk\gfile\elements\S.out"
# path = "D:\BaiduSyncdisk\gfile\elements\CO.out"
path="D:\BaiduSyncdisk\gfile\elements\H2.out"
path="D:\BaiduSyncdisk\gfile\elements\He2.out"

mol = Mol(reader=LogReader(path))

caler=density.Calculator(mol)

# nums=[]
# for i in range(5):
#     result=caler.atmDens(i+1)
#     print(i+1,np.sum(result))
#     nums.append(np.sum(result))
# print(sum(nums))

## 每个原子轨道电子密度加和应该为1
atmPos=caler.atmPos(1)
# print(atmPos.shape)

molPos,molWei=caler.molPos
# print(molPos.shape,molWei.shape)

# wfn=caler.wfnCaler.atoWfn(1,atmPos)
# den=wfn**2
# wei=caler.weights
# print(np.sum(den*wei))

# wfn=caler.wfnCaler.atoWfn(1,molPos)
# den=wfn**2
# wei=caler.a2mWeight(1)
# print(np.sum(den*wei)) # 相当于在原来正确的基础上添加了一些东西，肯定就不对了

## 分子轨道的电子密度
# print('-'*20)
# for obt in mol.O_obts:
#     wfn=caler.wfnCaler.obtWfn(obt,caler.points)
#     den=wfn**2
#     wei=caler.weights
#     print(np.sum(den*wei))

# print('-'*20)
# for obt in mol.O_obts:
#     wfn=caler.wfnCaler.obtWfn(obt,molPos)
#     den=wfn**2
#     wei=molWei
#     print(np.sum(den*wei))

# print('-'*20)
# qs=[]
# for atom in mol.atoms:
#     dens=caler.atmDens_ca(atom.idx)
#     q=np.sum(dens)
#     print(q)
#     qs.append(q)
# print(sum(qs))

# print('-'*20)
# qs=[]
# for atom in mol.atoms:
#     dens=caler.atmDens_cm(atom.idx)
#     q=np.sum(dens)
#     qs.append(q)
#     print(q)
# print(sum(qs))

coord=np.random.rand(20,3)*10
# print(coord)
# dens=caler.molDens_lib(coord)
# print(np.min(dens),np.max(dens))
# print(dens)
# print(mol.CM)

# obts=mol.O_obts
# nobt=len(obts)
# print(mol.CM[:,obts])
# print(mol.CM[:,:nobt])

dens=caler.molDens_lib(molPos)
print(np.sum(dens*molWei))

print('-'*40)
dens=caler.molDens_obt(molPos)
print(np.sum(dens*molWei))



# wfn=caler.wfnCaler.obtWfn(0,molPos)
# den=wfn**2
# wei=caler.a2mWeight(2)
# print(np.sum(den*wei))
