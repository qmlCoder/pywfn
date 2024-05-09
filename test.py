import numpy as np
from pywfn.base import Mol
from pywfn.reader import LogReader
from pywfn.atomprop import activity
from pywfn.maths import points_rotate
import pyvista as pv

# def normed(vect:np.ndarray):
#     length=np.linalg.norm()
#     return vect.copy()/length

# va=np.random.rand(3)
# vb=np.random.rand(3)
# va/=np.linalg.norm(va)
# vb/=np.linalg.norm(vb)

# cent=np.zeros(3)
# vc=np.cross(va,vb)
# vc/=np.linalg.norm(vc)
# print(np.linalg.norm(va))
# print(np.linalg.norm(vb))
# print(np.linalg.norm(vc))
# p=(vc+cent).reshape(1,3)
# axis=vb-va
# axis/=np.linalg.norm(axis)

# angles=np.linspace(0,np.pi,18,endpoint=True)
# points=np.zeros(shape=(len(angles),3))
# for a,ang in enumerate(angles):
#     res=points_rotate(p,cent,axis,ang).flatten()
#     points[a]=res

import time
t0=time.time()
name='NHC-wfn'
root=r"D:\BaiduSyncdisk\Articles\HFV\gfile\NHC"
molN=Mol(reader=LogReader(f"{root}/{name}-.log"))
mol0=Mol(reader=LogReader(f"{root}/{name}0.log"))
molP=Mol(reader=LogReader(f"{root}/{name}+.log"))

t1=time.time()
print(t1-t0)

t0=time.time()
caler=activity.Calculator()
caler.mols=[molN,mol0,molP]
result=caler.dirFukui(atms=[1])
t1=time.time()
print(t1-t0)

print('亲电性'.center(40,'-'))
for a,x,y,z,e,n in result:
    print(f'[{a:>3.0f},{x:>6.2f},{y:>6.2f},{z:>6.2f},{e:>6.2f}],')
print('亲核性'.center(40,'-'))
for a,x,y,z,e,n in result:
    print(f'[{a:>3.0f},{x:>6.2f},{y:>6.2f},{z:>6.2f},{n:>6.2f}],')

import linecache
path=r"D:\BaiduSyncdisk\Articles\HFV\gfile\M4\M4_wfn-.log"
line=linecache.getline(path,1)
print(line)
# start=5000000
# for i in range(start,start+20):
#     line=linecache.getline(path,i)
#     print(line)

# t0=time.time()
# lineNum=0
# with open(path,'r') as f:
#     for line in f:
#         lineNum+=1
# t1=time.time()
# print(t1-t0,lineNum)


# t0=time.time()
# lineNum=0
# with open(path,'r') as f:
#     while True:
#         line=f.readline()
#         if not line:break
#         lineNum+=1
# t1=time.time()
# print(t1-t0,lineNum)


# t0=time.time()
# with open(path,'r') as f:
#     for i  in range(lineNum):
#         line=linecache.getline(path,i)
# t1=time.time()
# print(t1-t0,lineNum)

# import re

# t0=time.time()
# string='123'
# for i in range(1000000):
#     '123' in string
# t1=time.time()
# print(t1-t0)

# t0=time.time()
# string='123'
# for i in range(1000000):
#     re.match('\d+',string)
# t1=time.time()
# print(t1-t0)