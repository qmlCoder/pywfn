import numpy as np

from pywfn import config
# from pywfn.base._mole import Mole
from pywfn.utils import printer

def cubeGrid(
    p0: np.ndarray, p1: np.ndarray, step: float, bord: float = 0
):
    """生成网格数据点,range:生成数据的范围
    getN:是否获取每个维度的数量
    getR:是否获取每个维度的长度
    """
    from pywfn import core
    assert isinstance(p0, np.ndarray), "必须是np.ndarray类型"
    p0 =p0.copy()- bord
    p1 =p1.copy()+ bord
    x0, y0, z0 = p0
    x1, y1, z1 = p1
    lx,ly,lz=dp=p1-p0
    assert dp.min() > 0 , "x输入的坐标范围错误"
    Nx=int((x1-x0)/step)+1
    Ny=int((y1-y0)/step)+1
    Nz=int((z1-z0)/step)+1
    xs=np.arange(Nx)
    ys=np.arange(Ny)
    zs=np.arange(Nz)
    XS,YS,ZS=np.meshgrid(xs,ys,zs)
    grids=np.array([
        XS.flatten(),
        YS.flatten(),
        ZS.flatten()
    ]).T
    grids=grids*step+p0
    return [Nx,Ny,Nz],grids

# 平面格点
def rectGrid(cent:np.ndarray,norm:np.ndarray,vx:np.ndarray,size:float)->tuple[list[int],np.ndarray]:
    """
    在空间中创建一个矩形区域,方便导出二维图像
    需要指定平面的法向量/z轴方向
    再指定平面的x轴方向，根据x和z轴计算y轴方向
    再根据z轴和y轴修复x轴方向
    """
    vz=norm.copy() 
    vz/=np.linalg.norm(vz) #平面的z轴朝向
    vx/=np.linalg.norm(vx) #平面的x轴朝向
    step=config.IMG_SPACE_STEP
    dxs=np.linspace(0,size,int(size/step))
    dys=np.linspace(0,size,int(size/step))
    
    nx=len(dxs)
    ny=len(dys)
    vy=np.cross(vx,vz)
    vx=np.cross(vz,vy)
    p0=cent-vx/2*size-vy/2*size
    # print(p0,cent,vx,vy)
    grid=[]
    for dx in dxs:
        for dy in dys:
            grid.append(p0+vx*dx+vy*dy)
    return [nx,ny],np.array(grid)

# 直线格点
def lineGrid(p0:np.ndarray,p1:np.ndarray,step:float):
    """获取一条线上的空间坐标"""
    vect=(p1-p0).astype(float)
    length=np.linalg.norm(vect).item()
    count=int(length/step)
    grid=[]
    vect/=length
    for i in np.linspace(0,length,count):
        grid.append(p0+i*vect)
    return np.array(grid)

def c2d_sph(xs,ys):
    """
    2D平面笛卡尔坐标转球谐坐标
    并不是所有的2D坐标都能转为球坐标
    """
    S=0.75
    ts=[]
    ps=[]
    for i,y in enumerate(ys):
        x=xs[i]
        p=y
        t=x/(np.cos(S*p))
        if t> np.pi:continue
        if t<-np.pi:continue
        ts.append(t)
        ps.append(p)
    ts=np.array(ts)
    ps=np.array(ps)
    return ts,ps

def sph_c3d(ts,ps,r=1):
    """球谐坐标转3D笛卡尔坐标"""
    xs=np.cos(ps)*np.cos(ts)*r
    ys=np.cos(ps)*np.sin(ts)*r
    zs=np.sin(ps)*r
    return xs,ys,zs

def sph_c2d(ts,ps):
    """球型坐标转2D坐标"""
    S=0.75
    xs=[]
    ys=[]
    for t,p in zip(ts,ps):
        x=t*np.cos(p*S)
        y=p
        xs.append(x)
        ys.append(y)
    xs=np.array(xs)
    ys=np.array(ys)
    return xs,ys

def CM2PM(CM:np.ndarray, obts: list[int], oe: int) -> np.ndarray:
    """
    根据系数矩阵构建密度矩阵
    CM:系数矩阵,如果是开壳层的话,列数是行数的两倍[n,n]/[n,2n]
    n:分子轨道占据电子数
    """
    assert type(obts)==list,"obts必须为列表"
    CMv=CM[:,obts].copy()
    CMh=np.transpose(CMv)
    PM = (CMv @ CMh) * oe
    return PM


def CM2PMs(CM, obts: list[int], oe: int):
    """
    构建三维密度矩阵，不要空轨道，形状为[占据轨道数,原子轨道数,原子轨道数]
    """
    A = (CM[:, obts].T)[:, :, np.newaxis]
    B = (CM[:, obts].T)[:, np.newaxis, :]
    return A @ B * oe  # 用矩阵乘法的形式直接构建矩阵可比逐元素计算多了


def vector_angle(a: np.ndarray, b: np.ndarray) -> float:  # 计算两向量之间的夹角
    """
    计算两向量之间的夹角
    a:o->a
    b:o->b
    """
    la = np.linalg.norm(a)
    lb = np.linalg.norm(b)
    assert la * lb > 0, "向量长度不应为0"

    cos = np.dot(a, b) / (la * lb)  # 应介于(-1,1)之间
    if cos > 1:
        cos = 1
    if cos < -1:
        cos = -1
    angle = np.arccos(cos) / np.pi
    return float(angle)


def linear_classify(points):  # 将向量分类转为表示角度的数值分类
    nv = points[-1]
    angles = np.array([vector_angle(each, nv) for each in points])
    y = np.abs(0.5 - angles)
    types = [[], []]
    idxs = [[], []]
    for i, each in enumerate(y):
        distance = np.abs(each - np.array([0, 0.5]))
        idx = np.argmin(distance)
        types[idx].append(each)
        idxs[idx].append(i)
    return idxs


def orbital_classify(atom, vectors, O_obts):
    """定义函数将轨道分类"""
    points = []
    orbitals = []
    keys = vectors.keys()
    for key in keys:
        atomidx = int(key.split("-")[0])
        orbitalidx = int(key.split("-")[1])
        if atomidx == atom and orbitalidx in O_obts and vectors[key] is not None:
            vector = vectors[key] / np.linalg.norm(vectors[key])
            points.append(vector)
            orbitals.append(orbitalidx)
    if len(points) > 0:
        points = np.array(points)
        idxs = linear_classify(points)
    else:
        idxs = [[], []]
    V, H = [[orbitals[id] for id in idx] for idx in idxs]
    return V, H


def get_aroundPoints(p, step):  #
    """
    计算空间中一个点周围六个点处的坐标
    p:[3,]坐标
    step:步长
    """
    arounds = np.array(
        [[-1, 0, 0], [+1, 0, 0], [0, -1, 0], [0, +1, 0], [0, 0, -1], [0, 0, +1]]
    )
    np.random.shuffle(arounds)
    return p + arounds * step

# def search_sp2(idx: int, mole: "Mole") -> int|None:
#     """
#     深度优先搜索方法寻找于指定原子相邻最近的sp2 C原子
#     """
#     atom = mole.atom(idx)

#     searchd = [idx]  # 已经搜索过的
#     searchs = [a for a in atom.neighbors]  # 将要搜索的

#     while len(searchs) > 0:
#         idx = searchs.pop(0)  # 弹出第一个
#         atom = mole.atom(idx)
#         if len(atom.neighbors) == 3:
#             return idx
#         searchn = [
#             a for a in mole.atom(idx).neighbors if a not in searchd
#         ]  # 新找到的原子的索引
#         searchs += searchn
#     printer.warn("没有找到sp2类型原子")
#     return None


def points_rotate(
    points: np.ndarray, center: np.ndarray, axis: np.ndarray, angle: float
):  # 旋转点
    """
    将一个点或一系列点绕空间中的某个轴旋转指定角度，输入和输出都是空间中点的绝对位置
    points:要旋转的点
    center:旋转中心
    axis:旋转轴
    angle:旋转角度
    """
    assert len(points.shape) == 2, "旋转点应该为二维数组"
    assert points.shape[1] == 3, "旋转点的形状应为[n,3]"
    axis /= np.linalg.norm(axis)
    nx,ny,nz=axis
    c=np.cos(angle)
    s=np.sin(angle)
    matrix = np.array([
        [nx**2*(1-c)+c,nx*ny*(1-c)+nz*s,nx*nz*(1-c)-ny*s],
        [nx*ny*(1-c)-nz*s,ny**2*(1-c)+c,ny*nz*(1-c)+nx*s],
        [nx*nz*(1-c)+ny*s,ny*nz*(1-c)-nx*s,nz**2*(1-c)+c],
    ])
    # print('旋转的点\n',points)
    # print('旋转矩阵\n',matrix)
    length0 = np.linalg.norm(points, axis=1)
    points = np.copy(points)
    points -= center
    points = (matrix@points.T).T
    points += center
    length1 = np.linalg.norm(points, axis=1)
    return points


def get_plane_by_3points(
    pos1: np.ndarray, pos2: np.ndarray, pos3: np.ndarray
) -> list[float]:
    """
    根据三个点计算平面
    平面方程：ax+by+cz+d=0

    """
    v21 = pos1 - pos2
    v23 = pos3 - pos2
    a, b, c = np.cross(v21, v23)
    x1, y1, z1 = pos1
    d = -a * x1 - b * y1 - c * z1
    return [a, b, c, d]


def get_distance_to_plane(params: list[float], pos: np.ndarray) -> float:
    """
    计算点到平面距离
    params：a,b,c,d
    """
    a, b, c, d = params
    x, y, z = pos
    distan = np.abs(a * x + b * y + c * z + d) / np.sqrt(a**2 + b**2 + c**2)
    return distan

def value2color(value:float,vmin:float,vmax:float,cmin:np.ndarray,cmax:np.ndarray):
    """将数值映射为颜色

    Args:
        value (float): 数值
        vmin (float): 数值下限
        vmax (float): 数值上限
        cmin (np.ndarray): 下限颜色
        cmax (np.ndarray): 上限颜色
    """
    k=(value-vmin)/(vmax-vmin)
    color=cmin+k*(cmax-cmin)
    return color

def gtf(grid,l,m,n)->np.ndarray:
    x=grid[:,0]
    y=grid[:,1]
    z=grid[:,2]
    r2=x**2+y**2+z**2
    alp=2.0
    facs=[1,1,3]
    fac=facs[l]*facs[m]*facs[n]
    ang=l+m+n
    Nm=(2*alp/np.pi)**(3/4)*np.sqrt((4*alp)**ang/fac)
    val=x**l * y**m * z**n * np.exp(-alp*r2)*Nm
    return val


def decomOrbitals(T:np.ndarray,coefs:np.ndarray,keeps:list):
    match len(coefs):
        case 1:
            return decomOrbitalS(T,coefs,keeps)
        case 3:
            return decomOrbitalP(T,coefs,keeps)
        case 6:
            return decomOrbitalD(T,coefs,keeps)
        case _:
            return coefs

def decomOrbitalS(T:np.ndarray,coefs:np.ndarray,keeps:list[int]):
    if keeps:
        return coefs
    else:
        return np.array([0.])

# 分解P轨道
def decomOrbitalP(T:np.ndarray,rcoefs:np.ndarray,keeps:list[int])->np.ndarray:
    """分解P轨道

    Args:
        T (np.ndarray): 基坐标，每一行代表一个方向
        rcoefs (np.ndarray): 原始函数空间的基函数系数
        keeps (list[int]): 保留的角动量

    Returns:
        np.ndarray: 分解之后的轨道系数
    """
    npos=3
    cords=np.random.rand(npos,3) #随机生成3个点[n,3],[3,3]
    wfn_1=np.zeros(shape=(npos,3)) #再这些点上的波函数数值
    wfn_1[:,0]=gtf(cords,1,0,0)
    wfn_1[:,1]=gtf(cords,0,1,0)
    wfn_1[:,2]=gtf(cords,0,0,1)

    wfn_2=np.zeros(shape=(npos,3))
    wfn_2[:,0]=gtf(cords@T,1,0,0)
    wfn_2[:,1]=gtf(cords@T,0,1,0)
    wfn_2[:,2]=gtf(cords@T,0,0,1)

    Mr=np.linalg.inv(wfn_2)@wfn_1
    Mi=np.linalg.inv(Mr)
    tcoefs=Mr@rcoefs # 根据函数空间基组1下的系数获取函数空间基组2下的系数
    for i in range(3):
        if i in keeps:continue
        tcoefs[i]=0.0
    fcoefs=Mi@tcoefs # 根据修改后的函数空间基组2下的系数得到函数空间基组1下的系数
    return fcoefs

# 分解D轨道
def decomOrbitalD(T:np.ndarray,rcoefs:np.ndarray,keeps:list[int]):
    npos=6
    cords=np.random.rand(npos,3)   # 随机生成6个点
    wfn_1=np.zeros(shape=(npos,6))
    wfn_1[:,0]=gtf(cords,2,0,0)
    wfn_1[:,1]=gtf(cords,0,2,0)
    wfn_1[:,2]=gtf(cords,0,0,2)
    wfn_1[:,3]=gtf(cords,1,1,0)
    wfn_1[:,4]=gtf(cords,1,0,1)
    wfn_1[:,5]=gtf(cords,0,1,1)

    wfn_2=np.zeros(shape=(npos,6))
    wfn_2[:,0]=gtf(cords@T,2,0,0)
    wfn_2[:,1]=gtf(cords@T,0,2,0)
    wfn_2[:,2]=gtf(cords@T,0,0,2)
    wfn_2[:,3]=gtf(cords@T,1,1,0)
    wfn_2[:,4]=gtf(cords@T,1,0,1)
    wfn_2[:,5]=gtf(cords@T,0,1,1)

    Mr=np.linalg.inv(wfn_2)@wfn_1
    Mi=np.linalg.inv(Mr)
    tcoefs=Mr@rcoefs
    for i in range(6):
        if i in keeps:continue
        tcoefs[i]=0.0
    fcoefs=Mi@tcoefs
    return fcoefs