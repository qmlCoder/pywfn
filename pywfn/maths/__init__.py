import numpy as np
import math
from pywfn import config
from pywfn import base
from pywfn.utils import printer
from pywfn.maths.gto import Gto


def gridPos(p0:np.ndarray,p1:np.ndarray,step:float,bord:float=0,getN:bool=False):
    """生成网格数据点,range:生成数据的范围
    getN:是否获取每个维度的数量
    getR:是否获取每个维度的长度
    """
    assert isinstance(p0,np.ndarray),'必须是np.ndarray类型'
    x0,y0,z0=p0-bord
    x1,y1,z1=p1+bord
    pos=[]
    xs=np.arange(x0,x1,step)
    ys=np.arange(y0,y1,step)
    zs=np.arange(z0,z1,step)
    Nx,Ny,Nz=len(xs),len(ys),len(zs)
    for x in xs:
        for y in ys:
            for z in zs:
                pos.append([x,y,z])
    
    if getN:
        return (Nx,Ny,Nz),np.array(pos,dtype=np.float32)
    else:
        return np.array(pos,dtype=np.float32)

def CM2PM(CM,obts:list[int],oe:int)->np.ndarray:
    """
    根据系数矩阵构建密度矩阵
    CM:系数矩阵,如果是开壳层的话,列数是行数的两倍[n,n]/[n,2n]
    n:分子轨道占据电子数
    """
    PMs=CM2PMs(CM,obts,oe)
    return np.sum(PMs,axis=0)

def CM2PMs(CM,obts:list[int],oe:int):
    """
    构建三维密度矩阵，不要空轨道，形状为[占据轨道数,原子轨道数,原子轨道数]
    """
    A=(CM[:,obts].T)[:,:,np.newaxis]
    B=(CM[:,obts].T)[:,np.newaxis,:]
    return A@B*oe #用矩阵乘法的形式直接构建矩阵可比逐元素计算多了

def vector_angle(a,b,trans=False): # 计算两向量之间的夹角
    '''计算两向量之间的夹角'''
    if np.linalg.norm(a)*np.linalg.norm(b)==0:
        raise
    value=np.dot(a,b)/(np.linalg.norm(a)*np.linalg.norm(b))
    if value>1:value=1
    if value<-1:value=-1
    angle=np.arccos(value)/math.pi
    if not trans:
        return angle
    else:
        return 0.5-abs(angle-0.5)
    
def get_normalVector(p1,p2,p3,linear=False):
    '''
    获取三点确定的平面的单位法向量
    根据两个向量也可以确定法向量
    '''
    vi=p3-p2
    vj=p1-p2
    n=np.cross(vi,vj) #法向量
    if np.linalg.norm(n)==0: # 此时说明三个原子在一条直线上，是标准的线型分子，所以不存在法向量
        return None
    if vector_angle(vi,vj,trans=True)<0.02 and linear:
        return None
    if vector_angle(n,config.BASE_VECTOR)>0.5:n*=-1
    return n/np.linalg.norm(n) # 返回单位向量

def linear_classify(points): #将向量分类转为表示角度的数值分类
    nv=points[-1]
    angles=np.array([vector_angle(each,nv) for each in points])
    y=np.abs(0.5-angles)
    types=[[],[]]
    idxs=[[],[]]
    for i,each in enumerate(y):
        distance=np.abs(each-np.array([0,0.5]))
        idx=np.argmin(distance)
        types[idx].append(each)
        idxs[idx].append(i)
    return idxs

def orbital_classify(atom,vectors,O_obts):
    '''定义函数将轨道分类'''
    points=[]
    orbitals=[]
    keys=vectors.keys()
    for key in keys:
        atomidx=int(key.split('-')[0])
        orbitalidx=int(key.split('-')[1])
        if atomidx==atom and orbitalidx in O_obts and vectors[key] is not None:
            vector=vectors[key]/np.linalg.norm(vectors[key])
            points.append(vector)
            orbitals.append(orbitalidx)
    if len(points)>0:
        points=np.array(points)
        idxs=linear_classify(points)
    else:
        idxs=[[],[]]
    V,H=[[orbitals[id] for id in idx] for idx in idxs]
    return V,H

def get_aroundPoints(p,step): # 
    '''
    计算空间中一个点周围六个点处的坐标
    p:[3,]坐标
    step:步长
    '''
    arounds=np.array([
        [-1,0,0],
        [+1,0,0],
        [0,-1,0],
        [0,+1,0],
        [0,0,-1],
        [0,0,+1]
    ])
    np.random.shuffle(arounds)
    return p+arounds*step

def get_extraValue(atom:"base.Atom",obt:int,valueType='max'):
    '''
    从指定位置开始,利用爬山算法寻找原子波函数极值
    maxPos:[3,]
    '''
    # p0=atom.coord.copy() # 起始点,p是原子坐标不变
    # p=p0.copy()
    p0=np.array([0.,0.,0.]).reshape(1,3)
    v0=atom.get_wfnv(p0,obt) # 计算原子坐标处的初始值
    step=0.1
    while True:
        aroundPs=get_aroundPoints(p0,step) # aroundPs:(n,3)
        aroundVs=atom.get_wfnv(aroundPs,obt)
        if valueType=='max' and np.max(aroundVs)>v0:
            maxID=np.argmax(aroundVs) #最大值的索引
            p0=aroundPs[maxID] # 最大值坐标
            v0=aroundVs[maxID] # 最大值
        elif valueType=='min' and np.min(aroundVs)<v0:
            minID=np.argmin(aroundVs)
            p0=aroundPs[minID]
            v0=aroundVs[minID]
        else:
            if step<=1e-6:
                return p0,v0 #
            else:
                step/=10

def search_sp2(idx:int,mol:"base.Mol")->int:
    """
    深度优先搜索方法寻找于指定原子相邻最近的sp2 C原子
    """
    atom=mol.atom(idx)
    
    searchd=[idx] #已经搜索过的
    searchs=[a.idx for a in atom.neighbors] #将要搜索的
    
    while len(searchs)>0:
        idx=searchs.pop(0) #弹出第一个
        atom=mol.atom(idx)
        if len(atom.neighbors)==3:return idx
        searchn=[a.idx for a in mol.atom(idx).neighbors if a.idx not in searchd] #新找到的原子的索引
        searchs+=searchn
    printer.warn('没有找到sp2类型原子')
    return None

def points_rotate(points:np.ndarray,axis_center:np.ndarray,axis_vector:np.ndarray,angle:float): # 旋转点
    """
    将一个点或一系列点绕空间中的某个轴旋转指定角度
    """
    axis_length=np.linalg.norm(axis_vector)
    axis_vector /=axis_length 
    a = np.cos(angle / 2.0)
    b, c, d = -axis_vector * np.sin(angle / 2.0)
    # 旋转矩阵
    matrix = np.array([
        [a**2 + b**2 - c**2 - d**2, 2 * (b*c - a*d), 2 * (b*d + a*c)],
        [2 * (b*c + a*d), a**2 - b**2 + c**2 - d**2, 2 * (c*d - a*b)],
        [2 * (b*d - a*c), 2 * (c*d + a*b), a**2 - b**2 - c**2 + d**2]
    ])
    points=np.copy(points)
    points-=axis_center
    points = np.dot(matrix, points.T).T
    points+=axis_center
    return points