"""
计算分子空间电子密度
整个分子只有网格点，但是没有权重
对于每个原子，根据自身的坐标及权重插值到整个空间的坐标及权重
？ 波函数值是否也要插值呢？

空间中一个格点的权重是否应该对所有分子都一致？
"""
from pywfn.base import Mol
from pywfn.spaceProp import wfnfunc
from pywfn.data import sphGrid
from functools import cached_property,lru_cache
from pywfn.spaceProp import dftgrid
from pywfn import maths
import numpy as np

class Calculator:
    def __init__(self,mol:Mol) -> None:
        self.mol=mol
        self.wfnCaler=wfnfunc.Calculator(mol)
        self.PM=self.mol.PM.copy()
        self.grids:np.ndarray=None
    
    def set_grid(self,grids:np.ndarray):
        self.grids=grids
        self.wfnCaler.set_grid(grids)
    
    # @lru_cache
    # def a2mWeight_(self,atm:int):
    #     """将原子格点的权重转为分子格点的权重""" # 径向和角度分别插值
    #     atmGrid,atmWeit=self.atmGrid(atm)
    #     atoms=self.mol.atoms
    #     natm=len(atoms)
    #     LM=self.mol.atoms.LM
    #     a2mPos=[]
    #     a2mWei=[]
    #     for p,gp in enumerate(atmGrid): # 对每个格点坐标进行循环
    #         S_u=np.ones(shape=(natm,natm))
    #         for i in range(natm):
    #             pi=atoms[i].coord
    #             ri=np.linalg.norm(gp-pi) # 格点与原子的距离
    #             for j in range(natm):
    #                 if i==j:continue # 相同原子，i=j，直接跳过了
    #                 pj=atoms[j].coord
    #                 rj=np.linalg.norm(gp-pj)
    #                 miu_ij=(ri-rj)/LM[j,i]
    #                 chi=atoms[i].radius/atoms[j].radius # 半径的比例，不受单位影响
                    
    #                 # Change μ(i,j) to ν(i,j)
    #                 if abs(chi-1)<1e-6: # 相同元素，半径相等
    #                     nu_ij=miu_ij
    #                 else:
    #                     u_ij=(chi-1)/(chi+1)
    #                     a_ij=u_ij/(u_ij**2-1)
    #                     if a_ij>0.5:a_ij=0.5
    #                     if a_ij<-0.5:a_ij=-0.5
    #                     nu_ij=miu_ij+a_ij*(1-miu_ij**2)

    #                 nu_ij=1.5*nu_ij-0.5*nu_ij**3
    #                 nu_ij=1.5*nu_ij-0.5*nu_ij**3
    #                 nu_ij=1.5*nu_ij-0.5*nu_ij**3
                    
    #                 S_u[j,i]=0.5*(1-nu_ij) # 和两原子的半径及间距有关
    #         wt=np.ones(natm)
    #         for i in range(natm):
    #             wt*=S_u[i,:]
    #         rat=wt[atm-1]/np.sum(wt) # 单个原子时，rat=1，相当于没做修改
    #         wi=atmWeit*rat
    #         a2mPos.append(gp)
    #         a2mWei.append(wi)
    #     return a2mPos,np.array(a2mWei)
    
    def molDens_obt(self,pos:np.ndarray):
        """根据分子轨道电子密度计算分子电子密度"""
        obts=self.mol.O_obts
        CM=self.mol.CM
        atms=self.mol.atoms.indexs
        dens=np.zeros(len(pos))
        for obt in obts:
            coefs=CM[:,obt]
            wfn=self.wfnCaler.obtWfn(obt,pos,atms,coefs)
            dens+=wfn**2*self.mol.oE
        return dens
    
    def molDens_lib(self):
        """使用Fortran库计算电子密度"""
        from pywfn.maths import flib
        ngrid=len(self.grids)
        nmat=self.mol.CM.shape[0]
        nobt=len(self.mol.O_obts)
        obts=self.mol.O_obts

        CM=self.mol.CM[:,obts].copy()

        wfns=self.wfnCaler.atoWfns
        dens=flib.molDens(ngrid,nmat,nobt,CM,wfns)
        return dens*self.mol.oE
        
    def atmDens(self,atm:int,grids:np.ndarray):
        """计算原子电子密度，使用分子空间坐标"""
        nmat=self.mol.CM.shape[0]
        atom=self.mol.atom(atm)
        dens=np.zeros(len(grids))
        u,l = atom.obtBorder
        self.wfnCaler.set_grid(grids)
        for i in range(u,l):
            wfn_i=self.wfnCaler.atoWfn(i)
            for j in range(nmat):
                wfn_j=self.wfnCaler.atoWfn(j)
                dens+=wfn_i*wfn_j*self.PM[i,j]
        return dens
    
    def rectValue(self,cent:np.ndarray,norm:np.ndarray,vx:np.ndarray,sx:float,sy:float):
        """生成图片文件"""
        nx,ny,grid=maths.rectGrid(cent,norm,vx,sx,sy)
        self.set_grid(grid)
        dens=self.molDens_lib()
        mat=dens.reshape(nx,ny)
        # plt.imshow(mat,cmap='Blues',vmin=0,vmax=1)
        # plt.savefig(f'density.jpg',bbox_inches='tight',pad_inches=0,dpi=300)
        return mat
    
    def lineValue(self,p0:np.ndarray,p1:np.ndarray):
        """获取一条直线上的数值"""
        grid=maths.lineGrid(p0,p1)
        self.set_grid(grid)
        dens=self.molDens_lib()
        return dens