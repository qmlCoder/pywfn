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
        return (self.nx,),self.grid

class RectGrid(Grid):
    def __init__(self):
        self.type:str='Rect'
        self.nx=1
        self.ny=1
        self.grid=np.array([[0,0,0]])
    
    def set(self,cent:np.ndarray,norm:np.ndarray,vx:np.ndarray,size:float):
        shape,grid=rectGrid(cent,norm,vx,size) # type: ignore
        self.grid=grid
        self.nx=shape[0]
        self.ny=shape[1]
        self.size=shape
        return self

    def get(self):
        return (self.nx,self.ny),self.grid

class CubeGrid(Grid):
    def __init__(self):
        self.type:str='Cube'
        self.nx=1
        self.ny=1
        self.nz=1
        self.grid=np.array([[0,0,0]])
        self.step=[1.,1.,1.]
    
    def set(self,p0:np.ndarray,p1:np.ndarray,step:float,bord=0):
        shape,grid=cubeGrid(p0,p1,step,bord)
        self.grid=grid
        self.nx=shape[0]
        self.ny=shape[1]
        self.nz=shape[2]
        self.size=shape
        self.step=[step,step,step]
        return self
    
    def get(self):
        return [self.nx,self.ny,self.nz],self.grid

class MapsGrid(Grid):
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
    
class VandGrid(Grid):
    def __init__(self):
        pass

    def set(self,mol):
        from pywfn.spaceprop import density
        densCaler=density.Calculator(mol)
        verts,faces=densCaler.vandSurf()
        self.grid=verts
        self.verts=verts
        self.faces=faces
        self.size=[len(verts),]
        return self

    def get(self):
        return self.grid

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
                plt.imshow(vals.reshape(*size),cmap='rwb')
                plt.savefig(path)
            case 'Maps':
                plt.imshow(vals.reshape(*size))
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
    
    def toJson(self,path:str,vals:np.ndarray,vand:VandGrid): ## 将范德华表面及对应的数值保存为json文件
        import json
        data={
            'verts':vand.verts.round(2).tolist(),
            'faces':vand.faces.round(2).tolist(),
            'value':vals.round(2).tolist()
        }
        Path(path).write_text(json.dumps(data,indent=4))
        
    def isoSurf(self,vals:np.ndarray,cube:CubeGrid,iosv:float):
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
        verts=np.concatenate(verts)
        faces=np.concatenate(faces)
        return verts,faces

from pywfn.shell import Shell

# 根据用户输入获取格点
def read_grid(shell:Shell,gidx:str,mol:Mol|None=None)->LineGrid|RectGrid|CubeGrid|VandGrid|MapsGrid:
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
            rect=RectGrid().set(cent,norm,vx,side)
            return rect
        case '3': # 空间采点
            p0  =shell.input.Float(tip='输入起点: ',count=3)
            p1  =shell.input.Float(tip='输入终点: ',count=3)
            step=shell.input.Float(tip='输入步长: ',count=1)[0] # type: ignore
            p0=np.array(p0)
            p1=np.array(p1)
            cube=CubeGrid().set(p0,p1,step)
            return cube
        case '4': # 分子格点
            assert mol is not None,"没有提供分子"
            p0,p1=mol.molBorder
            cube=CubeGrid().set(p0,p1,0.2,5)
            return cube
        case '5': # 地图
            from pywfn.spaceprop import MapsGrid
            r=shell.input.Float(tip='输入半径: ',count=1)[0] # type: ignore
            map=MapsGrid().set(r,100)
            return map
        case '6': # 范德华表面
            vand=VandGrid().set(mol)
            return vand
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
            obje:MapsGrid = obje # type: ignore
            vals=vals.flatten()
            vals[obje.idxs.flatten()] = 0.0
            caler.toImg(f'maps.png',vals,obje.size,'Maps')
        case 'Cube':
            assert obts is not None,"没有提供轨道索引"
            caler.toCub(f'cube.cub',vals,obje,obts) # type: ignore
        case 'Vand':
            caler.toJson(f'vand.json',vals,obje) # type: ignore

def onShell(shell:Shell):
    while True:
        printer.options('空间性质',{
            '1':'  波函数',
            '2':'电子密度',
            '3':'  静电势',
        })
        ctype=input('请选择计算的性质: ')
        mol=shell.input.Moles()[0]
        show_dict({'1':'直线采点','2':'平面采点','3':'空间采点','4':'分子空间','5':'地图映射','6':'范德华表面'})
        gidx=input('请选择采点方式: ')
        obje=read_grid(shell,gidx,mol)
        match ctype:
            case '1':
                from pywfn.spaceprop import wfnfunc
                caler=wfnfunc.Calculator(mol)
                atms=shell.input.Integ(tip='?输入要计算的原子: ')
                if not atms:atms=mol.atoms.atms
                obts=shell.input.Integ(tip='*输入要计算的轨道: ')
                wfns=caler.atmWfns(obje.grid,atms,obts)
                save_file(caler,obje,wfns,obts=obts)
            case '2':
                from pywfn.spaceprop import density
                caler=density.Calculator(mol)
                atms=shell.input.Integ(tip='?输入要计算的原子: ')
                if not atms:atms=mol.atoms.atms
                dens=caler.atmDens(obje.grid,atms)
                dens=dens.reshape(1,-1)
                save_file(caler,obje,dens)
            case '3':
                from pywfn.spaceprop import potential
                caler=potential.Calculator(mol)
                pots=caler.molPotential(obje.grid)
                save_file(caler,obje,pots)