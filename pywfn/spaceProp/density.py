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
from scipy.interpolate import griddata,LinearNDInterpolator
from functools import cached_property,lru_cache

import numpy as np

class Calculator:
    def __init__(self,mol:Mol) -> None:
        self.mol=mol
        
        self.points=sphGrid.gridData[:,:-1] # 网格的坐标
        self.weights=sphGrid.gridData[:,-1] # 网格的权重
        self.wfnCaler=wfnfunc.Calculator(mol)
        self.wfnCaler.molPos=self.molPos
    
    @cached_property
    def molPos(self):
        """整个分子的网格点坐标"""
        natm=len(self.mol.atoms)
        # 对格点畸形矫正
        DM=self.mol.atoms.DM
        print(DM)
        atoms=self.mol.atoms
        molPos=[]
        molWei=[]
        for a in range(natm):
            atmPos=self.atmPos(a+1)
            for p,gp in enumerate(atmPos): # 对每个格点坐标进行循环
                S_u=np.ones_like(DM)
                for i in range(natm):
                    pi=atoms[i].coord
                    r_i=np.linalg.norm(gp-pi) # 格点与原子的距离
                    
                    for j in range(natm):
                        if i==j:continue # 相同原子，i=j，直接跳过了
                        pj=atoms[j].coord
                        r_j=np.linalg.norm(gp-pj)

                        miu_ij=(r_i-r_j)/DM[j,i]
                        chi=atoms[i].radius/atoms[j].radius # 半径的比例，不受单位影响
                        
                        # Change μ(i,j) to ν(i,j)
                        if abs(chi-1)<1e-6: # 相同元素，半径相等
                            nu_ij=miu_ij
                        else:
                            uij=(chi-1)/(chi+1)
                            aij=uij/(uij**2-1)
                            if aij>0.5:aij=0.5
                            if aij<-0.5:aij=-0.5
                            nu_ij=miu_ij+aij*(1-miu_ij**2)

                        nu_ij=1.5*nu_ij-0.5*nu_ij**3
                        nu_ij=1.5*nu_ij-0.5*nu_ij**3
                        nu_ij=1.5*nu_ij-0.5*nu_ij**3
                        
                        S_u[j,i]=0.5*(1-nu_ij) # 和两原子的半径及间距有关
                wt=np.ones(natm)
                for i in range(natm):
                    wt*=S_u[i,:]
                rat=wt[a]/np.sum(wt) # 单个原子时，rat=1，相当于没做修改
                wi=self.weights[p]*rat
                molPos.append(gp)
                molWei.append(wi)
        return np.array(molPos),np.array(molWei)

    def atmPos(self,atm:int):
        """单个原子的网格点坐标"""
        return self.points.copy()+self.mol.atom(atm).coord
    
    @lru_cache
    def a2mWeight(self,atm:int): # 将原子的格点权重插值到分子中，但是权重是密度的权重而不是波函数的权重
        """单个原子在分子网格中的权重""" # 径向和角度分别插值
        atmPos=self.atmPos(atm)
        molPos=self.molPos.copy()-self.mol.atom(atm).coord.reshape(1,3) # 分子空间格点移动到以原子为中心
        radData=sphGrid.radData
        sphData=sphGrid.sphData
        # 根据插值获取权重
        rads=np.linalg.norm(molPos,axis=1) # 径向部分大小
        sphs=molPos/rads.reshape(-1,1) # 单位角度部分
        radw=np.interp(rads,radData[:,0],radData[:,1]) # 径向插值
        sphw=griddata(sphData[:,:-1],sphData[:,-1],sphs,method='linear',fill_value=0)
        a2mWeight=radw*sphw
        assert True not in np.isnan(a2mWeight),"权重计算不正确"
        return a2mWeight
    
    def a2mDens(self,atm:int):
        """将原子的电子密度插值到整个分子空间"""
        atmPos=self.atmPos(atm)
        molPos=self.molPos
        atmDens=self.atmDens(atm)
        a2mDens=griddata(atmPos,atmDens,molPos,method='linear')
        return a2mDens

    def molDens_obt(self,pos:np.ndarray):
        """计算根据分子轨道电子密度计算的分子电子密度"""
        obts = self.mol.O_obts
        dens=np.zeros(len(pos))
        for o in obts:
            wfn=self.wfnCaler.obtWfn(o,pos)
            dens+=wfn**2*self.mol.oE
        return dens
    
    def molDens_atm(self,pos:np.ndarray):
        """根据原子电子密度计算的分子密度"""
        dens=np.zeros(len(pos))
        for atom in self.mol.atoms:
            dens+=self.atmDens(atom.idx,pos)
        return dens
    
    def atmDens(self,atm:int):
        """计算原子电子密度，使用原子中心坐标"""
        atmPos=self.atmPos(atm)
        nmat=self.mol.CM.shape[0]
        atom=self.mol.atom(atm)
        dens=np.zeros(len(atmPos))
        atmDens=np.zeros(len(atmPos))
        u,l = atom.obtBorder
        for i in range(u,l):
            wfn_i=self.wfnCaler.atoWfn(i,atmPos)
            for j in range(nmat):
                wfn_j=self.wfnCaler.atoWfn(j,atmPos)
                dens=wfn_i*wfn_j*self.mol.PM[i,j]
                # assert True not in np.isnan(dens),"密度计算不正确"
                atmDens+=dens
        return atmDens*self.weights
    
    def atmDens2(self,atm:int):
        """计算原子电子密度，使用分子空间坐标"""
        molPos,molWei=self.molPos
        nmat=self.mol.CM.shape[0]
        atom=self.mol.atom(atm)
        dens=np.zeros(len(molPos))
        u,l = atom.obtBorder
        for i in range(u,l):
            wfn_i=self.wfnCaler.atoWfn(i,molPos)
            for j in range(nmat):
                wfn_j=self.wfnCaler.atoWfn(j,molPos)
                dens+=wfn_i*wfn_j*self.mol.PM[i,j]
                # assert True not in np.isnan(dens),"密度计算不正确"
        return dens*molWei
    
    def weiSqrt(self,weights:np.ndarray):
        """返回波函数的权重"""
        unit=weights/np.linalg.norm(weights) # 正负
        return np.sqrt(np.abs(weights))*unit