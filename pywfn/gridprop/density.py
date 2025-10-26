"""
计算分子空间电子密度
整个分子只有网格点，但是没有权重
对于每个原子，根据自身的坐标及权重插值到整个空间的坐标及权重
？ 波函数值是否也要插值呢？

空间中一个格点的权重是否应该对所有分子都一致？
"""
from pywfn.base.mole import Mole
from pywfn.gridprop import wfnfunc
from pywfn import gridprop
from pywfn.data import radDens
from pywfn import utils
import numpy as np
from pywfn import core

class Calculator(gridprop.SpaceCaler):
    def __init__(self,mole:Mole) -> None:
        self.mole=mole
        self.wfnCaler=wfnfunc.Calculator(mole)
        self.caler=core.gridprop.density.Calculator(mole.mole) # type: ignore # 核心计算器

    def mol_rho(self,grids:np.ndarray,level:int): # 以分子轨道的方式计算电子密度，可以计算每个分子轨道的电子密度
        """计算分子电子密度"""
        return self.caler.mol_rho(grids,level)
    
    def pi_rho(self,grids:np.ndarray):
        return self.caler.pi_rho(grids)
    
    def ato_rho(self,grids:np.ndarray,level:int):
        """计算所有原子轨道的电子密度"""
        return self.caler.ato_rho(grids,level)
    
    def fragDens(self,grids:np.ndarray,atms:list[int]|None=None):
        """计算指定原子的电子密度加和，使用分子空间坐标"""
        from pywfn import core
        if atms is None:atms=self.mole.atoms.atms
        dens=np.zeros(shape=(len(atms),len(grids)),dtype=float)
        xyzs,lmns,coes,alps=self.mole.basis.atoMap()
        atoDens=core.space.ato_rhos(grids,xyzs,lmns,coes,alps,self.mole.PM,0)[0] # type: ignore

        atoDens=np.array(atoDens)
        for a,atm in enumerate(atms):
            atom=self.mole.atom(atm)
            u,l = atom.obtBorder
            dens[a,:]=atoDens[:,u:l].sum(axis=1)
        return dens

    def molDens(self,grids:np.ndarray,level:int)->tuple[np.ndarray,np.ndarray,np.ndarray]:
        """使用Rust库计算电子密度

        Args:
            grids (np.ndarray): 要计算的格点
            level (int): 计算等级

        Returns:
            tuple[np.ndarray,np.ndarray,np.ndarray]: 电子密度、电子密度导数、电子密度Hessian
        """
        from pywfn import core
        CM=self.mole.coefs.get_CM('car')[:,self.mole.O_obts].copy() # 占据轨道的系数矩阵
        xyzs,lmns,coes,alps=self.mole.basis.atoMap()
        dens0,dens1,dens2=core.space.mol_rhos(grids,xyzs,lmns,coes,alps,CM,level) # type: ignore
        dens0=np.array(dens0)*self.mole.oE
        dens1=np.array(dens1)*self.mole.oE
        dens2=np.array(dens2)*self.mole.oE
        return dens0,dens1,dens2
    
    def obtDens(self,grids:np.ndarray,level:int,obt:int):
        """计算指定分子轨道的电子密度"""
        from pywfn import core
        oldCM=self.mole.coefs.get_CM('car')[:,self.mole.O_obts].copy() # 占据轨道的系数矩阵
        newCM=np.zeros_like(oldCM)
        newCM[:,obt]=oldCM[:,obt]
        xyzs,lmns,coes,alps=self.mole.basis.atoMap()

        dens0,dens1,dens2=core.space.mol_rhos(grids,xyzs,lmns,coes,alps,newCM,level) # type: ignore
        dens0=np.array(dens0)*self.mole.oE
        dens1=np.array(dens1)*self.mole.oE
        dens2=np.array(dens2)*self.mole.oE
        return dens0,dens1,dens2
    
    def piDens(self,grids:np.ndarray,level:int):
        from pywfn import core
        from pywfn.orbtprop import piobt
        obtCaler=piobt.Calculator(self.mole)
        CMp=obtCaler.project()
        xyzs,lmns,coes,alps=self.mole.basis.atoMap()
        dens0,dens1,dens2=core.space.mol_rhos(grids,xyzs,lmns,coes,alps,CMp,level) # type: ignore
        dens0=np.array(dens0)*self.mole.oE
        dens1=np.array(dens1)*self.mole.oE
        dens2=np.array(dens2)*self.mole.oE
        return dens0,dens1,dens2
    
    def proDens(self,grids:np.ndarray): # 第一种方式，使用插值计算
        """计算前体分子电子密度 promol"""
        dens=np.zeros((len(grids)))
        for atom in self.mole.atoms:
            radius=np.linalg.norm(grids-atom.coord,axis=1)
            dens+=radDens.get_radDens(atom.atomic,radius)
        return dens
    
    def proDens_v2(self,grid:np.ndarray,level=0): #第二种版本，使用公式计算
        """计算前体分子电子密度及梯度 promol"""
        dens0=np.zeros((len(grid),))
        dens1=np.zeros((len(grid),3))
        dens2=np.zeros((len(grid),3,3))
        for atom in self.mole.atoms:
            rho0,rho1,rho2=radDens.get_radDens_v2(atom.atomic,grid,level)
            dens0+=rho0
            dens1+=rho1
            dens2+=rho2
        return dens0,dens1,dens2
            
    def RDG(self,grids:np.ndarray,pro:bool=False): # reduced density gradient
        assert utils.chkArray(grids,[None,3]),"格点的形状应为: (n,3)"
        if pro: # 如果使用预分子
            dens0,dens1,_=self.proDens_v2(grids,1) # type: ignore
        else: # 使用真实电子密度及梯度
            dens0,dens1,_=self.molDens(grids,2)
        L=np.linalg.norm(dens1,axis=1)
        K=1/(2*(3*np.pi**2)**(1/3))
        rdg=K*L/dens0**(4/3)
        return rdg
    
    def IRI(self,grids:np.ndarray,pro:bool=False,alp=1.1): # 
        assert utils.chkArray(grids,[None,3]),"格点的形状应为: (n,3)"
        if pro: # 如果使用预分子
            dens0,dens1,_=self.proDens_v2(grids,1) # type: ignore
        else:
            dens0,dens1,_=self.molDens(grids,1)
        L=np.linalg.norm(dens1,axis=1)
        # return L/dens0**a
        
        vals=L/dens0**alp
        # print(np.min(vals),np.max(vals))
        return vals
    
    def signl2rho(self,grids:np.ndarray,pro:bool=False):
        assert utils.chkArray(grids,[None,3]),"格点的形状应为: (n,3)"
        if pro: # 如果使用预分子
            dens0,dens1,dens2=self.proMolDens_v2(grids,1) # type: ignore
        else:
            dens0,dens1,dens2=self.molDens(grids,2)
        resVals=np.zeros_like(dens0)
        for i,each in enumerate(dens2):
            eigVals=np.linalg.eigvalsh(each)
            eigVal=eigVals[1]
            if eigVal==0:
                resVals[i]=0
                continue
            resVals[i]=eigVal/np.abs(eigVal)*dens0[i]
            if resVals[i]<-0.05:resVals[i]=-0.05
            if resVals[i]> 0.05:resVals[i]= 0.05
            
        return resVals
    
    def dletaDens(self,grids:np.ndarray):
        """差分电子密度"""
        proDens=self.proDens(grids)
        molDens=self.molDens(grids,0)[0]
        return molDens-proDens

    def laplace(self,grids:np.ndarray,pro:bool=False):
        """计算电子密度的拉普拉斯"""
        assert utils.chkArray(grids,[None,3]),"格点的形状应为: (n,3)"
        if pro: # 如果使用预分子
            dens0,dens1,dens2=self.proMolDens_v2(grids,1) # type: ignore
        else:
            dens0,dens1,dens2=self.molDens(grids,2)
        return dens2[:,0,0]+dens2[:,1,1]+dens2[:,2,2]
    
    def piMolDens(self,grid:np.ndarray):
        """计算分子空间π电子密度"""
        pass

    def vandSurf(self):
        """返回分子的范德华表面"""
        p0,p1=self.mole.spaceBorder
        shape,grids=gridprop.CubeGrid().set_v1(p0,p1,0.3,4).get()
        vals=self.molDens(grids,0)[0]
        verts,faces,_,_=self.isoSurf(shape,grids,vals,0.04)
        assert verts is not None,"未找到范德华表面"
        assert faces is not None,"未找到范德华表面"
        return verts,faces
    
    def AIM(self):
        """AIM分析，寻找键临界点"""
        def show_type(Hessian):
            eigval,eigvec=np.linalg.eigh(Hessian)
            rank=[1 if abs(e)>1e-6 else  0 for e in eigval]
            syms=[1 if e>0         else -1 for e in eigval]
            # print(sum(rank),sum(syms))
            return sum(rank),sum(syms)
        crits=[]
        for bond in self.mole.bonds:
            a1,a2=bond.ats
            pos0=(self.mole.atom(a1).coord+self.mole.atom(a2).coord)/2.0
            # print(pos0)
            for i in range(100):
                dens0,dens1,dens2=self.molDens(pos0.reshape(1,3),level=2) #计算电子密度，电子密度梯度和电子密度Hessian矩阵
                Hessian=dens2[0] #取出Hess矩阵
                # eigvals,eigvecs=np.linalg.eigh(Hessian) #计算Hess矩阵的特征值和特征向量
                if np.abs(np.linalg.det(Hessian)) < 1e-10:
                    x,y,z=pos0.tolist() # type: ignore
                    r,s=show_type(Hessian)
                    crits.append([x,y,z,r,s])
                    break
                if np.linalg.norm(dens1[0]) < 1e-6:
                    # print(dens1[0])
                    r,s=show_type(Hessian)
                    x,y,z=pos0.tolist() # type: ignore
                    crits.append([x,y,z,r,s])
                    break
                delta = np.linalg.solve(Hessian, dens1[0])
                pos0 -= delta
                # print(i,pos0,np.linalg.norm(dens1[0]))
        return crits