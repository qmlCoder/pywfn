"""
计算原子的各种方向
- 原子法向量
- p轨道方向
- 原子周围一圈方向
"""
from pywfn.base.mole import Mole
from pywfn import maths
from pywfn.maths import vector_angle,points_rotate,cubeGrid
from pywfn import config
import re
import numpy as np
from functools import lru_cache
from collections import deque
from pywfn import core

maxWeaves={} # 记录已经计算过的原子最大值方向

class Calculator:
    def __init__(self,mole:Mole) -> None:
        self.mole=mole
        self.noNorms=[] # 已经计算过的原子法向量
        self.normals={} # 存储原子的法向量
        self.caler=core.atomprop.direction.Calculator(mole.mole) # type: ignore # 核心计算器

    def normal_vector(self,atm:int):
        """计算原子的法向量"""
        return self.caler.normal_vector(atm)
    
    def get_normal(self,atm:int):
        """计算原子的法向量"""
        return self.caler.get_normal(atm)
    
    def local_coord_system(self,atm:int,neb:int|None=None):
        """计算原子的局部坐标系"""
        return self.caler.local_coord_system(atm,neb)
    
    def reactions(self,atm:int):
        """计算原子可能的反应方向"""
        return self.caler.reactions(atm)
    
    @lru_cache
    def maxWeave(self,atm:int,obt:int,sym:str)->np.ndarray|None:
        """
        计算原子波函数最大值方向
        atm:原子
        obt:分子轨道
        sym:轨道的符号,正则表达式
        """
        def get_vals(pos,step): # 计算周围一圈的波函数值
            delPos=np.array([
                [+step,0,0],[-step,0,0],
                [0,+step,0],[0,-step,0],
                [0,0,+step],[0,0,-step]
            ])
            vals=np.zeros(6) # 只计算六个点的波函数值
            for i in idxs: # 遍历所有符合要求的轨道
                coef=CM[i,obt]
                wfn=wfnCaler.obt_wfn(pos+delPos,obt,0)[0]
                vals+=coef*wfn # 分子轨道=组合系数*原子轨道

            maxIdx=np.argmax(vals) # 最大值所在的索引
            maxVal=np.max(vals) # 最大值
            maxDir=delPos[maxIdx]
            return maxVal,maxDir
        from pywfn.gridprop import wfnfunc
        wfnCaler=wfnfunc.Calculator(self.mole) # 波函数计算器
        step=0.1

        CM=self.mole.CM.copy()
        nmat=CM.shape[0]
        idxs=[] # 符合要求的轨道索引
        for i in range(nmat):
            if self.mole.atoAtms[i]!=atm:continue #不是指定的原子不算
            obtSym=self.mole.atoSyms[i]
            if re.match(sym,obtSym) is None:continue # 符号不对不算
            idxs.append(i)


        pos0=self.mole.atom(atm).coord.copy()
        # print(pos0)
        val0=-np.inf #起初为负无穷
        for i in range(1000):
            maxVal,maxDir=get_vals(pos0,step)
            if maxVal>val0 and maxVal>0:
                val0= maxVal
                pos0+=maxDir
            elif step>1e-10:
                step/=10
            else:
                break
        
        direction=pos0-self.mole.atom(atm).coord
        length=np.linalg.norm(direction)
        if length<1e-8:
            return None
        else:
            return direction/length
    
    
    def sphAround(self)->np.ndarray|None: # 直接用dft的角度网格算了
        """计算周围一圈的球形范围"""
        from pywfn.data import lebedev
        lebedev.LD0074()
        
    def bases(self)->dict[int,np.ndarray]:
        """
        计算每个原子的局部基坐标[natm,3,3]，每一列代表一个基向量
        """
        atms=self.mole.heavyAtoms
        # base=np.zeros(shape=(natm,3,3))
        base={}
        for i,atm in enumerate(atms):
            atom=self.mole.atom(atm)
            nebs=atom.neighbors # 相邻的原子
            norm=self.normal_vector(atom.idx) # 原子法向量
            if norm is None:
                if atom.is_linear:
                    atm1,atm2=atom.neighbors
                elif len(atom.neighbors)==1:
                    atm1=atom.neighbors[0]
                    atm2=atm
                else:
                    raise ValueError(f"原子{atm}无法给定法向量")
                vx=self.mole.atom(atm1).xyz-self.mole.atom(atm2).xyz
                vx/=np.linalg.norm(vx)
                randVect=np.random.rand(3)
                randVect/=np.linalg.norm(randVect)
                norm=np.cross(vx,randVect)
                norm/=np.linalg.norm(norm)
                self.normals[atm]=norm
            # print(atm,norm)
            cent=self.mole.xyzs[[e-1 for e in nebs],:].copy().mean(axis=0) # 邻居原子坐标的平均值
            vect=atom.coord-cent # 原子到邻居原子的向量
            
            vz=norm/np.linalg.norm(norm) # 单位法向量
            # 决定vx
            if atom.is_linear:
                # print(f'原子{atom.idx}是线性的')
                idx1,idx2=nebs
                vx=self.mole.atom(idx1).coord-self.mole.atom(idx2).coord
                vx/=np.linalg.norm(vx)
            else:
                dist=np.linalg.norm(vect)
                # print(f'原子{atom.idx}到其邻心距离{dist}',vect)
                if dist>1e-1: # 两个点不重合
                    vx=vect/dist
                    if vector_angle(vz,vx)>0.5:vz*=-1
                else:
                    anyv=np.random.rand(3) # 任意向量
                    vx=np.cross(vz,anyv)
                    vx/=np.linalg.norm(vx)
            vy=np.cross(vz,vx)
            vy/=np.linalg.norm(vy)
            base[atm]=np.array([vx,vy,vz]).T # 原子局部坐标系的基坐标向量
        return base
    
    def hmoBases(self)->dict[int,np.ndarray]:
        """休克尔分子轨道法原子的方向"""
        from pywfn.orbtprop import popul
        NM=popul.Calculator(self.mole).piEleMatDecom('atom',False)
        eles=np.sum(NM,axis=0)
        piIdx=0 # 第一个pi分子轨道的索引
        for i,each in enumerate(eles):
            if abs(each-2)>0.5:continue
            piIdx=i
            break

        atms=self.mole.heavyAtoms
        base={}
        for i,atm in enumerate(atms):
            atom=self.mole.atom(atm)
            norm=self.normal_vector(atom.idx) # 原子法向量
            if norm is None:continue
            u,l=atom.obtBorder
            coefs=self.mole.CM[u:l,piIdx]
            syms=self.mole.atoSyms[u:l]
            piCoefs=[]
            for s,c in zip(syms,coefs):
                if 'P' not in s:continue
                piCoefs.append(c)
            piDir=np.array(piCoefs).reshape(-1,3).sum(axis=0)
            if vector_angle(piDir,norm)>0.5:norm*=-1 # type: ignore
            anyv=np.random.rand(3) # 任意向量
            vz=norm/np.linalg.norm(norm) # 单位法向量
            vx=np.cross(vz,anyv)
            vx/=np.linalg.norm(vx)
            vy=np.cross(vz,vx)
            vy/=np.linalg.norm(vy)
            base[atm]=np.array([vx,vy,vz]).T # 原子局部坐标系的基坐标向量
        return base

    # 线性原子环轴一圈的方向
    def linearAtomCircle(self,atm:int):
        from pywfn.maths import points_rotate
        atom=self.mole.atom(atm)
        if not atom.is_linear:return None
        a1,a2=atom.neighbors
        bondDir=self.mole.atom(a1).coord-self.mole.atom(a2).coord
        bondDir/=np.linalg.norm(bondDir)
        anyDir=np.random.rand(3)
        anyDir/=np.linalg.norm(anyDir)
        vertDir=np.cross(bondDir,anyDir)
        vertDir/=np.linalg.norm(vertDir)
        points=(atom.coord+vertDir).reshape(1,3)
        center=atom.coord.copy().reshape(1,3)
        dirs=[]
        for i in range(36):
            rotated_point=points_rotate(points,center,bondDir,np.pi*2/36*i)
            dirs.append(rotated_point.flatten()-atom.coord)
        return np.array(dirs)