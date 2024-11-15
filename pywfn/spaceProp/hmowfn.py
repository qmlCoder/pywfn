"""
休克尔分子轨道法轨道的可视化

将空间格点转换为原子局部坐标系的点，计算这些点上的波函数值
"""

from pywfn.base import Mol
from pywfn.reader import LogReader
from pywfn.atomprop import direction

from pywfn.maths.mol import hmo

import numpy as np

def hmoWfn(Z:int,grid:np.ndarray):
    a0=0.589
    n=1/(4*np.sqrt(2*np.pi))*(Z/a0)**(5/2)
    r=np.linalg.norm(grid,axis=1) # 到原子核的距离
    wfn=n*np.exp(-Z*r*0.5/a0)*grid[:,-1]
    return wfn


class Calculator:
    def __init__(self,mol:Mol):
        self.mol=mol
        self.grid=np.zeros(shape=(1,3))
        BM,es,CM,occs=hmo(self.mol)
        self.CM=CM

    def set_grid(self,grid:np.ndarray):
        self.grid=grid
    
    def atom_base(self):
        """计算每个原子的局部基坐标[natm,3,3]，每一列代表一个基向量
        """
        dirCaler=direction.Calculator(self.mol)
        atms=self.mol.heavyAtoms
        natm=len(atms)
        # base=np.zeros(shape=(natm,3,3))
        base={}
        for i,atm in enumerate(atms):
            atom=self.mol.atom(atm)
            norm=dirCaler.normal(atom.idx) # 原子法向量
            anyv=np.random.rand(3) # 任意向量
            vz=norm/np.linalg.norm(norm) # 单位法向量
            vx=np.cross(vz,anyv)
            vx/=np.linalg.norm(vx)
            vy=np.cross(vz,vx)
            vy/=np.linalg.norm(vy)
            base[atm]=np.array([vx,vy,vz]).T # 原子局部坐标系的基坐标向量
            # print(atm,np.linalg.norm(base[i],axis=0))
        return base

    def ato_wfns(self):
        """计算所有的原子轨道波函数
        """
        base=self.atom_base()
        atms=self.mol.heavyAtoms
        natm=len(atms)
        wfns=np.zeros(shape=(natm,len(self.grid)))
        for i,atm in enumerate(atms):
            atom=self.mol.atom(atm)
            grid=(self.grid-atom.coord)@base[atm] # 将空间格点转为以原子为中心的坐标
            wfn=hmoWfn(atom.atomic,grid)
            wfns[i]=wfn
        return wfns
    
    def ato_wfn(self,atm:int):
        base=self.atom_base()
        # atms=self.mol.heavyAtoms
        # natm=len(atms)
        # wfns=np.zeros(shape=(natm,len(self.grid)))
        atom=self.mol.atom(atm)
        grid=(self.grid-atom.coord)@base[atm] # 将空间格点转为以原子为中心的坐标
        wfn=hmoWfn(atom.atomic,grid)
        return wfn

    
    def obtWfn(self,obt:int)->np.ndarray:
        wfns=self.ato_wfns()
        atms=self.mol.heavyAtoms
        natm=len(atms)
        wfn=np.zeros(len(self.grid))
        for i,atm in enumerate(atms):
            coe=self.CM[i,obt]
            wfn+=coe*wfns[i]
        return wfn
    
    def sph_2DWfn(self,r:float,obt:int,size:int)->np.ndarray:
        """
        计算球坐标映射为2d坐标的波函数
        返回二维矩阵
        """
        # R=2.29/0.529177
        
        R=r
        print('半径:',R)
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
        idx=np.argwhere(np.logical_or(lons>np.pi,lons<-np.pi)) # 不符合条件的角度找出来
        

        # 计算波函数值
        self.set_grid(grids)
        wfns=self.obtWfn(obt)
        wfns[idx]=0
        matrix=wfns.reshape(ny,nx)
        return matrix

    def toCub(self,path:str,obts:list[int],pos0:np.ndarray,shape:list[int],step:list[float]):
        """导出成.cub文件

        Args:
            path (str): 导出的文件位置
            obts (list[int]): 渲染的轨道
            pos0 (np.ndarray): 格点起始位置
            shape (list[int]): 格点的形状
            step (list[float]): 格点的步长
        """
        from pywfn.writer import CubWriter
        vals=[self.obtWfn(obt) for obt in obts]
        syms=self.mol.atoms.syms
        xyzs=self.mol.coords
        nums=shape
        writer=CubWriter(syms,xyzs,[o+1 for o in obts],pos0,nums,step,vals)
        writer.save(path)