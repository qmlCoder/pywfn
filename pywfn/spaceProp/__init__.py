"""
定义一个平面的方式有很多种，但是本质上对应的格点都是一样的，因此可以用统一的类来实现
"""

import numpy as np
from numpy.typing import NDArray
from pywfn.maths import lineGrid,rectGrid,cubeGrid
from pywfn.utils import printer
from pywfn.maths import march
from pywfn.base import Mol
from pathlib import Path
from pywfn.maths import flib,rlib
import time

class Grid:
    def __init__(self):
        self.type:str=''
        self.size:list[int]  # 格点形状
        self.grid:np.ndarray

class LineGrid(Grid):
    def __init__(self):
        self.type:str='Line'
        self.nx=1
        self.grid=np.array([[0,0,0]])
    
    def set(self,p0:np.ndarray,p1:np.ndarray,step:float):
        grid=lineGrid(p0,p1,step)
        self.grid=grid
        self.nx=len(grid)
        self.size=[self.nx,]
        return self
    
    def get(self):
        return [self.nx,],self.grid

class RectGrid(Grid):
    def __init__(self):
        self.type:str='Rect'
        self.grids=np.array([[0,0,0]])
    
    def set_v1(self,cent:np.ndarray,norm:np.ndarray,vx:np.ndarray,size:float):
        shape,grids=rectGrid(cent,norm,vx,size) # type: ignore
        self.grids=grids
        self.shape=shape
        return self

    def set_v2(self,p0:np.ndarray,p1:np.ndarray,p2:np.ndarray): # 根据三个原子的坐标生成网格
        vx=p1-p0
        vn=p1-p2
        vx/=np.linalg.norm(vx)
        vn/=np.linalg.norm(vn)
        norm=np.cross(vx,vn)
        norm/=np.linalg.norm(norm)
        vy=np.cross(vx,norm)
        dist=np.linalg.norm(p1-p0).item() #两点之间的距离
        print('dist',dist)
        print('vx',vx)
        print('vy',vy)
        print('norm',norm)
        # start=p0-vx*4-vy*4 #设置起始点位置
        cent=p0+vx*dist/2+vy*dist/2
        size=dist+8
        print(cent,norm,vx,size)
        shape,grids=rectGrid(cent,norm,vx,size)
        self.shape=shape
        self.grids=grids
        return self

    def get(self):
        return self.shape,self.grids

class CubeGrid(Grid):
    def __init__(self):
        self.type:str='Cube'

        self.grid=np.array([[0,0,0]])
        self.step=[1.,1.,1.]
    
    def set_v1(self,p0:np.ndarray,p1:np.ndarray,step:float,bord:float=0.0):
        shape,grid=cubeGrid(p0,p1,step,bord)
        self.grid=grid

        self.shape=shape
        self.step=[step,step,step]
        return self
    
    def get(self):
        return self.shape,self.grid

class EarthGrid(Grid):
    def __init__(self):
        self.type:str='Maps'
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
        self.size=[ny,nx]
        return self

    def get(self)->tuple[list[int],np.ndarray]:
        return self.size,self.grid

class SpaceCaler:
    def __init__(self):
        self.mol:Mol
    
    def toImg(self,path:str,vals:np.ndarray,size:list[int],_type:str): #传入来的数据是一维的
        import matplotlib.pyplot as plt
        match _type:
            case 'Line':
                xs=np.arange(vals.shape[1])
                for each in vals:
                    plt.plot(xs,each)
                plt.savefig(path)
            case 'Rect':
                vmin=vals.min()
                vmax=vals.max()
                plt.imshow(vals.reshape(*size),cmap='bwr',vmin=vmin,vmax=vmax)
                plt.savefig(path)
            case 'Maps':
                vmin=vals.min()
                vmax=vals.max()
                plt.imshow(vals.reshape(*size),cmap='bwr',vmin=vmin,vmax=vmax)
                plt.savefig(path)

    def toCub(self,path:str,vals:np.ndarray,cube:CubeGrid,obts:list[int]): ## 保存为cub文件，cub文件是不是只能用来存储波函数？
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
        writer=CubWriter()
        writer.syms=syms
        writer.xyzs=xyzs
        writer.obts=obts
        writer.pos0=pos0
        writer.size=size
        writer.step=step
        writer.vals=vals

        writer.save(path)
    
    def toObj(self,path:str,verts:np.ndarray,faces:np.ndarray): ## 保存为obj文件
        from pywfn.writer import ObjWriter
        writer=ObjWriter()
        writer.verts=verts
        writer.faces=faces
        writer.save(path)
    
    @staticmethod
    def isoSurf(shape:list[int],grids:np.ndarray,vals:np.ndarray,isov:float,limit:tuple[float,float]|None=None,gt:bool=True):
        # faces=[]
        # voxelData  =march.grids2voxel(shape,grids,vals)
        # verts=march.voxel2verts(voxelData,isov,limit,gt) # 顶点坐标，每三个点代表一个面，包含很多重复的点
        print('提取等值面',limit,gt)
        # verts:np.ndarray=flib.marchCube(grids,vals,shape,isov)[1] # type: ignore 
        verts=rlib.march_cube(shape,grids.tolist(),vals.tolist(),isov) # type: ignore
        verts=np.array(verts)
        # print(verts)
        if verts is not None:
            print('vert.shape',verts.shape)
            # verts,faces=flib.vertsMerge(verts,0.1) # 合并顶点
            # faces=faces.reshape(-1,3)
            faces=np.arange(len(verts)).reshape(-1,3)
            # if not gt:
            #     faces[:,0],faces[:,2]=faces[:,2],faces[:,0]
            return verts,faces
        
        else:
            return None,None

from pywfn.shell import Shell

# 根据用户输入获取格点
def read_grid(shell:Shell,gidx:str,mol:Mol|None=None)->LineGrid|RectGrid|CubeGrid|EarthGrid:
    match gidx:
        case '1': # 直线采点
            p0  =shell.input.Float(tip='输入起点: ',count=3)
            p1  =shell.input.Float(tip='输入终点: ',count=3)
            step=shell.input.Float(tip='输入步长: ',count=1)[0] # type: ignore
            p0=np.array(p0)
            p1=np.array(p1)
            line=LineGrid().set(p0,p1,step)
            return line
        case '2': # 平面采点
            cent=shell.input.Float(tip='输入中心: ',count=3)
            norm=shell.input.Float(tip='输入法线: ',count=3)
            vx  =shell.input.Float(tip='输入 x轴: ',count=3)
            side=shell.input.Float(tip='输入边长: ',count=1)[0] # type: ignore
            cent=np.array(cent)
            norm=np.array(norm)
            vx=np.array(vx)
            rect=RectGrid().set_v1(cent,norm,vx,side)
            return rect
        case '3': # 空间采点
            p0  =shell.input.Float(tip='输入起点: ',count=3)
            p1  =shell.input.Float(tip='输入终点: ',count=3)
            step=shell.input.Float(tip='输入步长: ',count=1)[0] # type: ignore
            p0=np.array(p0)
            p1=np.array(p1)
            cube=CubeGrid().set_v1(p0,p1,step)
            return cube
        case '4': # 分子格点
            assert mol is not None,"没有提供分子"
            p0,p1=mol.molBorder
            cube=CubeGrid().set_v1(p0,p1,0.2,5)
            return cube
        case '5': # 地图
            from pywfn.spaceprop import EarthGrid
            r=shell.input.Float(tip='输入半径: ',count=1)[0] # type: ignore
            map=EarthGrid().set(r,100)
            return map
        case _:
            assert False,"无效的格点类型"


def show_dict(dict:dict):
    for k,v in dict.items():
        print(f'{k}: {v}; ',end='')
    print('\n')

def save_file(caler:SpaceCaler,obje:Grid,vals:np.ndarray, obts:list[int]|None=None): # type: ignore
    match obje.type:
        case 'Line':
            caler.toImg(f'line.png',vals,obje.size,'Line')
        case 'Rect':
            caler.toImg(f'rect.png',vals,obje.size,'Rect')
        case 'Maps':
            obje:EarthGrid = obje # type: ignore
            vals=vals.flatten()
            vals[obje.idxs.flatten()] = 0.0
            caler.toImg(f'maps.png',vals,obje.size,'Maps')
        case 'Cube':
            assert obts is not None,"没有提供轨道索引"
            caler.toCub(f'cube.cub',vals,obje,obts) # type: ignore