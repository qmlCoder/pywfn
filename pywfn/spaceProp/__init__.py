"""
定义一个平面的方式有很多种，但是本质上对应的格点都是一样的，因此可以用统一的类来实现
"""

import numpy as np
from numpy.typing import NDArray
from pywfn.maths import lineGrid,rectGrid,cubeGrid
from pywfn.utils import printer
from pywfn.maths import march
from pywfn.base import Mol

class Line:
    def __init__(self):
        self.nx=1
        self.grid=np.array([[0,0,0]])
    
    def set(self,p0:np.ndarray,p1:np.ndarray,step:float):
        grid=lineGrid(p0,p1,step)
        self.grid=grid
        self.nx=len(grid)
        return self
    
    def get(self):
        return (self.nx,),self.grid

class Rect:
    def __init__(self):
        self.nx=1
        self.ny=1
        self.grid=np.array([[0,0,0]])
    
    def set(self,cent:np.ndarray,norm:np.ndarray,vx:np.ndarray,size:float):
        size,grid=rectGrid(cent,norm,vx,size)
        self.grid=grid
        self.nx=size[0]
        self.ny=size[1]
        return self

    def get(self):
        return (self.nx,self.ny),self.grid

class Cube:
    def __init__(self):
        self.nx=1
        self.ny=1
        self.nz=1
        self.grid=np.array([[0,0,0]])
        self.step=1.0
    
    def set(self,p0:np.ndarray,p1:np.ndarray,step:float):
        size,grid=cubeGrid(p0,p1,step)
        self.grid=grid
        self.nx=size[0]
        self.ny=size[1]
        self.nz=size[2]
        self.size=size
        self.step=step
        return self
    
    def get(self):
        return (self.nx,self.ny,self.nz),self.grid

class Map:
    def __init__(self):
        self.nx=1
        self.ny=1
        self.grid=np.array([[0,0,0]])
        self.idxs=np.array([[0]])

    def set(self,r:float,size:int):
        R=r
        nx=size*2
        ny=size
        # 地图坐标
        xr=np.linspace(-2*R*2**0.5,2*R*2**0.5,nx)
        yr=np.linspace(-R*2**0.5,R*2**0.5,ny)
        xs,ys=np.meshgrid(xr,yr)
        xs=xs.flatten()
        ys=ys.flatten()
        # 经纬度
        theta=np.arcsin(ys/(R*2**0.5))
        lons=np.pi*xs/(2*R*2**0.5*np.cos(theta)) # 经度
        lots=np.arcsin((2*theta+np.sin(2*theta))/np.pi) # 纬度

        # 三维空间坐标
        grids=np.zeros(shape=(nx*ny,3))
        grids[:,1]=R*np.sin(lots) # y轴作为南北极

        grids[:,0]=R*np.cos(lots)*np.cos(lons)
        grids[:,2]=R*np.cos(lots)*np.sin(lons)
        idxs=np.argwhere(np.logical_or(lons>np.pi,lons<-np.pi)) # 不符合条件的角度找出来
        self.nx=nx
        self.ny=ny
        self.grid=grids
        self.idxs=idxs
        return self

    def get(self):
        return (self.nx,self.ny),self.grid,self.idxs
    

class SpaceCaler:
    def __init__(self):
        self.mol:Mol
    
    def toImg(self,path:str,vals:np.ndarray):
        import matplotlib.pyplot as plt
        ndim=len(vals.shape) #数据的维度
        if ndim==1:
            plt.plot(vals)
            plt.savefig(path)
        elif ndim==2:
            plt.imshow(vals)
            plt.savefig(path)

    def toCub(self,path:str,vals:np.ndarray,cube:Cube,obts:list[int]): ## 保存为cub文件，cub文件是不是只能用来存储波函数？
        """导出成.cub文件

        Args:
            path (str): 导出的文件位置
            vals (np.ndarray): 导出的数据，三维的数据
            obts (list[int]): 渲染的轨道
        """
        from pywfn.writer import CubWriter
        size,grid=cube.get()
        syms=self.mol.atoms.syms
        xyzs=self.mol.coords
        pos0=grid[0]
        step=cube.step
        writer=CubWriter(syms,xyzs,obts,pos0,size,step*3,vals)
        writer.save(path)
    
    def toObj(self,path:str,vals:np.ndarray,cube:Cube,iosv:float): ## 保存为obj文件
        from pywfn.writer import ObjWriter
        verts=[]
        faces=[]
        cubeData=march.cord2cube(cube.size,cube.grid,vals)
        vertP,vertN=march.cube2vert(cubeData,iosv) # 顶点坐标，每三个点代表一个面，包含很多重复的点
        npos=0
        for vert in [vertP,vertN]:
            if vert is None:continue
            ## 不过滤
            # verts.append(vert)
            # face=np.arange(1,len(vert)+1).reshape(-1,3)
            # faces.append(face)
            ## 过滤
            vert,face=march.filtVerts(vert) # 过滤掉重复的点
            face+=npos
            verts.append(vert)
            faces.append(face)
            npos+=len(vert)
        # faceP=np.array(len(vertP)//3,3)
        # faceN=np.array(len(vertN)//3,3)
        verts=np.concatenate(verts)
        faces=np.concatenate(faces)
        writer=ObjWriter(verts,faces)
        writer.save(path)


from pywfn.shell import Shell

def get_grid(shell:Shell,vtype:str):
    match vtype:
        case '直线':
            from pywfn.spaceProp import Line
            p0  =shell.input.Float(tip='输入起点: ',count=3)
            p1  =shell.input.Float(tip='输入终点: ',count=3)
            step=shell.input.Float(tip='输入步长: ',count=1)
            p0=np.array(p0)
            p1=np.array(p1)
            line=Line().set(p0,p1,step)
            size,grid=line.get()
            return size,grid,None,line
        case '平面':
            from pywfn.spaceProp import Rect
            cent=shell.input.Float(tip='输入中心: ',count=3)
            norm=shell.input.Float(tip='输入法线: ',count=3)
            vx  =shell.input.Float(tip='输入x轴: ',count=3)
            side=shell.input.Float(tip='输入边长: ',count=1)[0]
            cent=np.array(cent)
            norm=np.array(norm)
            vx=np.array(vx)
            rect=Rect().set(cent,norm,vx,side)
            size,grid=rect.get()
            return size,grid,None,rect
        case '空间':
            from pywfn.spaceProp import Cube
            p0  =shell.input.Float(tip='输入起点: ',count=3)
            p1  =shell.input.Float(tip='输入终点: ',count=3)
            step=shell.input.Float(tip='输入步长: ',count=1)
            p0=np.array(p0)
            p1=np.array(p1)
            cube=Cube().set(p0,p1,step)
            size,grid=cube.get()
            return size,grid,None,cube

def get_fileType(shell:Shell):    
    types={'1':'png','2':'cub','3':'obj'}
    show_dict(types)
    idx=shell.input.Integ(tip='请选择文件的类型: ')[0]
    return types[f'{idx}']

def save_file(caler:SpaceCaler,ftype:str,vals:np.ndarray,path:str,
              cube:Cube=None,obts:list[int]=None,isov:float=None):
    match ftype:
        case 'img':
            caler.toImg(path,vals)
        case 'cub':
            caler.toCub(path,vals,cube,obts)
        case 'obj':
            caler.toObj(path,vals,cube,isov)

def show_dict(dict:dict):
    for k,v in dict.items():
        print(f'{k}: {v}; ',end='')
    print('\n')
def onShell(shell:Shell):
    while True:
        printer.options('空间性质',{
            '1':'波函数',
            '2':'电子密度',
            '3':'静电势',
        })
        ctype=input('请选择计算的性质: ')
        vdict={'1':'直线','2':'平面','3':'空间','4':'地图'}
        show_dict(vdict)
        vtype=input('请选择数据的形式: ')
        size,grid,idxs,obje=get_grid(shell,vdict[vtype])
        match ctype:
            case '1':
                from pywfn.spaceProp import wfnfunc
                mol=shell.input.Moles()[0]
                caler=wfnfunc.Calculator(mol)
                atms=shell.input.Integ(tip='?输入要计算的原子: ')
                if not atms:atms=mol.atoms.atms
                obts=shell.input.Integ(tip='*输入要计算的轨道: ')
                wfns=caler.atmWfns(grid,atms,obts)
                ftype=get_fileType(shell)
                save_file(caler,ftype,wfns,f'波函数_{vdict[vtype]}.{ftype}',cube=obje,obts=obts)
            case '2':
                from pywfn.spaceProp import density
                mol=shell.input.Moles()[0]
                caler=density.Calculator(mol)
                atms=shell.input.Integ(tip='?输入要计算的原子: ')
                dens=caler.atmDens(grid,atms)
                ftype=get_fileType(shell)
                save_file(caler,ftype,dens,f'电子密度_{vdict[vtype]}.{ftype}',cube=obje,obts=[0])
                