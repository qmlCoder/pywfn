"""
实现march cube算法
根据三维坐标生成体素数据
输入体素数据，输出三维坐标
等值面提取，可以提取大于某个值或者小于某个值的等值面，大于或小于决定了朝向
设置函数阈值，不满足阈值的直接跳过，提高效率
python做等值面提取还是太慢了，得想办法用fortran来实现
"""
from pywfn.data.march import marchData
from pywfn import utils

from pathlib import Path
import json
import numpy as np

bonds = [[0, 1], [1, 2], [2, 3], [0, 3], [1, 5], [2, 6], [3, 7], [0, 4], [4, 5], [5, 6], [6, 7], [7, 4]]; ## 八个顶点确定的键
fdata:dict[str,list[list[int]]]=marchData

def grids2voxel(shape:list[int],grids:np.ndarray,values:np.ndarray)->np.ndarray: # 格点数据转体素数据
    """将坐标和值转换为体素数据 [nx,ny,nz,4] (x,y,z,v)，方便获取周围点坐标

    Args:
        shape (tuple[int,int,int]): 体素的形状
        cords (np.ndarray): 所有体素的坐标
        values (np.ndarray): 体素对应的数值

    Returns:
        np.ndarray: 体素数据
    """
    # print(values)
    nx,ny,nz=shape
    voxel=np.zeros(shape=(nx,ny,nz,4))
    n=0
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                voxel[i,j,k,:3]=grids[n]
                voxel[i,j,k,3]=values[n]
                # print(i,j,k,values[n])
                n+=1
    return voxel

def voxel2verts(cube:np.ndarray,isov:float,limit:tuple[float,float]|None=None,gt:bool=True)->np.ndarray|None: # 根据体素数据生成等值面的顶点 
    """获取体素数据的等值面坐标

    Args:
        cube (np.ndarray): 体素数据
        isov (float): 等值面数值

    Returns:
        np.ndarray: 顶点坐标
    """
    nx,ny,nz,_=cube.shape
    verts=[]
    # verticesN=[]
    for i in range(nx-1):
        for j in range(ny-1):
            for k in range(nz-1):
                coords,values=get_around(cube,i,j,k)
                # print(values)
                minVal=values.min()
                maxVal=values.max()
                if limit is not None:
                    if maxVal<limit[0]:continue
                    if minVal>limit[1]:continue
                # print(i,j,k,minVal,maxVal)
                if minVal< isov <maxVal: # 在正方形格点的极值之间代表在正方形内
                    keyP=''.join([str(int(v >  isov)) for v in values])
                    # print(keyP,int(keyP,base=2)+1,len(marchData[keyP]))
                    posP=get_vertices(keyP,values,coords,isov,gt)
                    verts.append(posP)
    if verts:
        verts=np.concatenate(verts)
    else:
        verts=None

    return verts

def filtVerts(verts:np.ndarray): # 过滤掉位置相同的顶点
    """顶点去重/顶带你合并"""
    assert utils.chkArray(verts,[None,3]),"verts must be a 2d array with shape (n,3)"
    # 计算每两个顶点之间的距离矩阵
    distMat = np.linalg.norm(verts[:, None, :] - verts[None, :, :], axis=-1)
    oldFaces=np.arange(len(verts)) # 原本的面索引
    # 找到距离小于1e-5的顶点对
    

    onMap={1:1} #新旧顶点之间的映射
    nvert=len(verts)
    filts=[verts[0]] #过滤出的点
    index=2
    for i in range(nvert):
        if i==0:continue
        nfilt=len(filts) #当前过滤出的点的数量
        start=0 if nfilt<100 else nfilt-100
        points=np.array(filts[start:]) #已经确定过的点
        npoint=verts[i,:].reshape(-1,3)
        dists=np.linalg.norm(points-npoint,axis=1) #当前点与之前点的距离
        minval=np.min(dists)
        minidx=np.argmin(dists)
        if minval<1e-2 and minidx!=len(filts)-1: # 如果距离小于1e-5，且不是当前点，则说明是重复点
            onMap[i+1]=start+minidx.item()+1 # 将当前点映射到重复点
        else:
            onMap[i+1]=index
            filts.append(verts[i])
            index+=1
    faces=[]
    for i in range(1,nvert+1,3):
        faces.append([onMap[i],onMap[i+1],onMap[i+2]])
    filts=np.array(filts)
    faces=np.array(faces)
    # print(f'过滤前:{len(verts)}个点，过滤后:{len(filts)}个点,index={index}')
    return filts,faces

def reguVerts(verts:np.ndarray,faces:np.ndarray): # 将每个顶点向周围顶点平均位置移动
    pass


def get_around(cube,i,j,k):
    didx = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]] # 周围8个点的相对索引
    values=[]
    coords=[]
    # print(i,j,k)
    for n,idx in enumerate(didx):
        di,dj,dk=idx
        i_ = i+di
        j_ = j+dj
        k_ = k+dk

        x=cube[i_,j_,k_,0]
        y=cube[i_,j_,k_,1]
        z=cube[i_,j_,k_,2]
        v=cube[i_,j_,k_,3]
        values.append(v)
        coords.append([x,y,z])
        # print(di,dj,dk)
        # print(f'{n:>3}{x:>10.4f}{y:>10.4f}{z:>10.4f}{v:>10.4f}{i_:>3}{j_:>3}{k_:>3}')
    values=np.array(values)
    coords=np.array(coords)
    return coords,values

def get_vertices(key:str,values:np.ndarray,coords:np.ndarray,isov:float,gt:bool):
    faces=fdata[key]
    verts=[]
    for f,face in enumerate(faces):
        # print(f,face)
        for p in range(3): #每一个面对应的三个点
            e=face[p]
            a,b=bonds[e]
            pa,pb=coords[a],coords[b]
            va,vb=values[a],values[b]
            # print(f,p,e,a,b,pa,pb,va,vb)
            if va==vb:
                k=0.5
            else:
                k = (isov - va)/(vb - va)
            dv = pb-pa
            point = pa+k*dv
            verts.append(point)
            # print(len(verts),k,point,va,vb)
    verts=np.array(verts)
    return verts