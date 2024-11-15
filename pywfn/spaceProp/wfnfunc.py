"""
计算分子轨道空间波函数
基函数的波函数
原子轨道的波函数
原子的波函数
"""
from pywfn.base import Mol
from pywfn import maths
import numpy as np
from pywfn.maths import flib
from pywfn.maths import cubeGrid
from pywfn.spaceProp import lutils
from pywfn.maths import march
Array=np.ndarray

class Calculator:
    def __init__(self,mol:Mol) -> None:
        self.mol=mol
        self.molPos=np.zeros((1,3)) # 初始坐标设为原点
        self.grids:np.ndarray=None # 格点数据
        self.wfns:np.ndarray=None
        self.CM=self.mol.CM.copy()
        self.atms=self.mol.atoms.atms
    
    def set_grid(self,grids:np.ndarray):
        assert len(grids.shape)==2,"grids必须为二维"
        assert grids.shape[1]==3,"grids必须为3列"
        self.grids=grids
        self.wfns=None # 重新设置格点的时候计算出的波函数要清零
    
    def obtWfn(self,obt:int):
        """
        计算分子轨道的波函数，为原子轨道的线性组合
        obt：分子轨道指标
        """
        assert self.grids is not None,"没有设置格点数据"
        coefs=self.CM[:,obt] # 轨道系数
        idxs=[]
        for atm in self.atms:
            u,l=self.mol.atom(atm).obtBorder
            idxs+=list(range(u,l))
        wfn=np.zeros(self.grids.shape[0])
        for c,coef in enumerate(coefs):
            if c not in idxs:continue
            wfn+=coef*self.atoWfn(c)
        return wfn
    
    def atoWfn(self,i:int):
        """
        计算原子轨道的波函数，形成在轨道线性组合之前
        pos,分子的空间坐标
        i:原子轨道指标
        """
        wfn=self.atoWfns[i,:]
        return wfn
    
    @property
    def atoWfns(self):
        """计算所有原子轨道"""
        assert self.grids is not None,"没有设置格点数据"
        if self.wfns is None:
            ngrid=self.grids.shape[0]
            nmat=self.mol.CM.shape[0]
            # cords=self.mol.coords.copy()
            atms = self.mol.obtAtms
            shls = self.mol.obtShls
            syms = self.mol.obtSyms
            angs = self.mol.obtAngs
            basDict = self.mol.basis.dict
            coords=[]
            expl=[]
            coel=[]
            ncgs=[]
            cmax=0
            for i in range(nmat):
                atm=atms[i]
                atom=self.mol.atom(atm)
                atomic=atom.atomic
                coords.append(atom.coord)
                shl=shls[i]
                ang=angs[i]
                key=f'{atomic}-{shl}-{ang}'
                dat=np.array(basDict[key])
                exps=dat[:,0]
                coes=dat[:,1]
                expl.append(exps)
                coel.append(coes)
                if len(coes)>cmax:cmax=len(coes)
                ncgs.append(len(coes))
            
            coords=np.array(coords)
            expa=np.zeros(shape=(nmat,cmax))
            coea=np.zeros(shape=(nmat,cmax))
            for i,ncg in enumerate(ncgs):
                expa[i,:ncg]=expl[i]
                coea[i,:ncg]=coel[i]

            lmns = [self.mol.basis.sym2lmn(sym) for sym in syms]
            lmns=np.array(lmns)
            ncgs=np.array(ncgs)
            
            self.wfns=flib.cgfs(ngrid,self.grids,nmat,coords,cmax,ncgs,expa,coea,lmns)
        return self.wfns
        
    def rectValue(self,cent:np.ndarray,norm:np.ndarray,vx:np.ndarray,sx:float,sy:float,atms:list[int],obt:int):
        """生成图片文件"""
        nx,ny,pos=maths.rectGrid(cent,norm,vx,sx,sy)
        wfns=np.zeros(shape=len(pos))
        nbas=self.mol.CM.shape[0]
        obtAtms=self.mol.obtAtms
        for i in range(nbas):
            if obtAtms[i] not in atms:continue
            coef=self.mol.CM[i,obt]
            wfns+=self.atoWfn(i,pos)*coef
        mat=wfns.reshape(nx,ny)
        return mat
    
    def sph_2DWfn(self,r:float,obt:int,size:int)->np.ndarray:
        """
        计算球坐标映射为2d坐标的波函数
        返回二维矩阵
        """
        # R=2.29/0.529177
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
        idx=np.argwhere(np.logical_or(lons>np.pi,lons<-np.pi)) # 不符合条件的角度找出来
        

        # 计算波函数值
        self.set_grid(grids)
        wfns=self.obtWfn(obt)
        wfns[idx]=0
        matrix=wfns.reshape(ny,nx)
        return matrix
    
    def lineValue(self,p0:np.ndarray,p1:np.ndarray,atms:list[int],obt:int):
        """获取一条直线上的数值"""
        grid=maths.lineGrid(p0,p1)
        wfns=np.zeros(shape=len(grid))
        nbas=self.mol.CM.shape[0]
        obtAtms=self.mol.obtAtms
        for i in range(nbas):
            if obtAtms[i] not in atms:continue
            coef=self.mol.CM[i,obt]
            wfn=self.atoWfn(i,grid)*coef
            wfns+=wfn
        return wfns
    
    def rectValue(self): # 在一个矩形内的波函数值
        grid=maths.rectGrid()
        pass

    def cubeValue(self,obt:int): # 在一个正方体内的波函数值
        p0,p1=lutils.get_molBorder(self.mol) # 获取分子边界
        shape,grid=cubeGrid(p0,p1,0.1,1)# 生成网格坐标
        self.set_grid(grid)
        wfns=self.obtWfn(obt).reshape(*shape)
        return wfns
    
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
    
    def toObj(self,path:str,wfns:np.ndarray): ## 保存为obj文件
        from pywfn.writer import ObjWriter
        vals=wfns.flatten()
        shape=wfns.shape
        cube=march.cord2cube(shape,self.grids,vals)
        vertP,vertN=march.cube2vert(cube,0.01)
        filtP,faceP=march.filtVerts(vertP)
        filtN,faceN=march.filtVerts(vertN)
        verts=np.concatenate([filtP,filtN])
        faces=np.concatenate(faceP,faceN)
        writer=ObjWriter(verts,faces)
        writer.save(path)


def meshgrid(xr,yr):
    xs,ys=np.meshgrid(xr,yr)
    xs=xs.flatten()
    ys=ys.flatten()
    return xs,ys