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

from pywfn.data import sphGrid
weight = sphGrid.gridData[:, -1]
coords = sphGrid.gridData[:, :3]

class Calculator:
    def __init__(self,mol) -> None:
        self.mol:Mol=mol
        self.molPos=np.zeros((1,3)) # 初始坐标设为原点
        self.grids:np.ndarray=None
        self.wfns:np.ndarray=None
        self.CM=self.mol.CM.copy()
        self.atms=self.mol.atoms.atms
    
    def set_grid(self,grids:np.ndarray):
        self.grids=grids
        self.wfns=None # 重新设置格点的时候计算出的波函数要清零
    
    def obtWfn(self,obt:int):
        """
        计算分子轨道的波函数，为原子轨道的线性组合
        obt：分子轨道指标
        coefs：线性组合系数
        atms：可以自定义原子
        CM：可自定义系数矩阵
        """
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
        # print('atoWfn',pos[0,:])
        # atms = self.mol.obtAtms
        # shls = self.mol.obtShls
        # syms = self.mol.obtSyms
        # lmns = [self.mol.basis.sym2lmn(sym) for sym in syms]
        
        # lmn = lmns[i]
        # atm = atms[i]
        # shl = shls[i]
        # ang = sum(lmn)
        # atmic = self.mol.atom(atm).atomic
        # basis = self.mol.basis.get(atmic, shl, ang)
        # exps = [b.exp for b in basis]
        # coes = [b.coe for b in basis]
        # coord=self.mol.atom(atm).coord
        # # R2 = np.sum(pos_**2, axis=1)
        # wfn = self.mol.gto.cgf(exps, coes, lmn, grids, coord)  # 空间坐标-以原子为中心的坐标

        # wfns=self.atoWfns(grids)
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
    
    def sph_2DWfn(self,r:float,obt:int)->np.ndarray:
        """
        计算球坐标映射为2d坐标的波函数
        返回二维矩阵
        """
        S=0.75
        nx,ny=1000,500

        xr=np.linspace(-np.pi,np.pi,nx)    # x取值范围
        yr=np.linspace(-np.pi/2,np.pi/2,ny)# y取值范围
        xs,ys=meshgrid(xr,yr)
        # 转换为球坐标
        ps=ys
        ts=xs/(np.cos(S*ys))
        # 记下不符合角度范围的索引
        idx=np.argwhere(np.logical_or(ts>np.pi,ts<-np.pi))
        
        # 生成3D笛卡尔坐标
        grids=np.zeros(shape=(nx*ny,3))
        grids[:,1]=np.cos(ps)*np.cos(ts)*r 
        grids[:,0]=np.cos(ps)*np.sin(ts)*r
        grids[:,2]=np.sin(ps)*r
        
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

def meshgrid(xr,yr):
    xs,ys=np.meshgrid(xr,yr)
    xs=xs.flatten()
    ys=ys.flatten()
    return xs,ys