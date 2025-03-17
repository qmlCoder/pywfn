from pywfn.base import Mole

import numpy as np

def get_molBorder(mol:Mole):
    coords=mol.coords
    p0=coords.min(axis=0)
    p1=coords.max(axis=0)
    p0-=0.1
    p1+=0.1
    return p0,p1

def sphArea(r:float):
    cords=[]
    cords.append([0,0,1])
    for i in np.linspace(0,np.pi,20)[1:-1]: #仰角
        for j in np.linspace(0,2*np.pi,40)[:-1]: #水平角
            z=np.cos(i)
            x=np.sin(i)*np.cos(j)
            y=np.sin(i)*np.sin(j)
            cords.append([x,y,z])
    cords.append([0,0,-1])
    return np.array(cords)*r

def vdeFace(mol:Mole): # 计算范德华表面
    
    cords=[]
    for i,atom in enumerate(mol.atoms):
        sphps=sphArea(atom.radius)
        print(atom,atom.radius)
        for eac in sphps+atom.coord:
            dists=np.linalg.norm(mol.coords-eac,axis=1) # 原子周围每个点到所有原子的距离
            inside=True
            for j,atomj in enumerate(mol.atoms):
                dist=dists[j]
                if dist<atomj.radius:
                    inside=False
            if inside:cords.append(eac)
    return np.array(cords)
