"""
实现march cube算法
根据三维坐标生成体素数据
输入体素数据，输出三维坐标
"""
from pywfn.config import ROOT_DATA
from pywfn.data.march import marchData
from pathlib import Path
import json
import numpy as np

bonds = [[0, 1], [1, 2], [2, 3], [0, 3], [1, 5], [2, 6], [3, 7], [0, 4], [4, 5], [5, 6], [6, 7], [7, 4]]; ## 八个顶点确定的键
data:dict[str,list[list[int]]]=marchData

def cord2cube(shape:tuple[int,int,int],cords:np.ndarray,values:np.ndarray)->np.ndarray:
    """将坐标和值转换为体素数据 [nx,ny,nz,4] (x,y,z,v)，方便获取周围点坐标

    Args:
        shape (tuple[int,int,int]): 体素的形状
        cords (np.ndarray): 所有体素的坐标
        values (np.ndarray): 体素对应的数值

    Returns:
        np.ndarray: 体素数据
    """
    nx,ny,nz=shape
    cube=np.zeros(shape=(nx,ny,nz,4))
    n=0
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                cube[i,j,k,:3]=cords[n]
                cube[i,j,k,3]=values[n]
                n+=1
    return cube

def cube2vert(cube:np.ndarray,isov:float)->tuple[np.ndarray|None,np.ndarray|None]: # 根据体素数据生成等值面的顶点
    """获取体素数据的等值面坐标

    Args:
        cube (np.ndarray): 体素数据
        isov (float): 等值面数值

    Returns:
        np.ndarray: 顶点坐标
    """
    nx,ny,nz,_=cube.shape
    verticesP=[]
    verticesN=[]
    for i in range(nx-1):
        for j in range(ny-1):
            for k in range(nz-1):
                coords,values=get_around(cube,i,j,k)
                minVal=values.min()
                maxVal=values.max()
                # print(i,j,k,minVal,maxVal)
                if minVal< isov<maxVal:
                    keyP=''.join([str(int(v>isov)) for v in values])
                    posP=get_vertices(keyP,values,coords,isov)
                    verticesP.append(posP)

                if minVal<-isov<maxVal:
                    keyN=''.join([str(int(v<-isov)) for v in values])
                    posN=get_vertices(keyN,values,coords,-isov)
                    verticesN.append(posN)
    if verticesP:
        verticesP=np.concatenate(verticesP)
    else:
        verticesP=None
    if verticesN:
        verticesN=np.concatenate(verticesN)
    else:
        verticesN=None
    return verticesP,verticesN

def filtVerts(verts:np.ndarray):
    """删除重复的顶点"""
    
    print('verts',verts.shape)
    print(verts[0])
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
        if minval<1e-3 and minidx!=len(filts)-1: # 如果距离小于1e-5，且不是当前点，则说明是重复点
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
    # print(filts)
    # print(f'过滤前:{len(verts)}个点，过滤后:{len(filts)}个点,index={index}')
    return filts,faces


def get_around(cube,i,j,k):
    didx = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]] # 周围8个点的相对索引
    values=[]
    coords=[]
    for idx in didx:
        di,dj,dk=idx
        i_=i+di
        j_=j+dj
        k_=k+dk

        x=cube[i_,j_,k_,0]
        y=cube[i_,j_,k_,1]
        z=cube[i_,j_,k_,2]
        v=cube[i_,j_,k_,3]
        values.append(v)
        coords.append([x,y,z])
    values=np.array(values)
    coords=np.array(coords)
    return coords,values

def get_vertices(key:str,values:np.ndarray,coords:np.ndarray,isov:float):
    faces=data[key]
    verts=[]
    for f,face in enumerate(faces):
        for p in range(3): #每一个面对应的三个点
            e=face[p]
            a,b=bonds[e]
            pa,pb=coords[a],coords[b]
            va,vb=values[a],values[b]
            k=(isov-va)/(vb-va)
            dv=pb-pa
            point=pa+k*dv
            verts.append(point)
    verts=np.array(verts)
    return verts