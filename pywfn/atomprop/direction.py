"""
计算原子的各种方向
- 原子法向量
- p轨道方向
- 原子周围一圈方向
"""
from pywfn.base import Mol
from pywfn import maths
from pywfn.maths import vector_angle,points_rotate,cubeGrid
from pywfn import config
import re
import numpy as np
from functools import lru_cache

maxWeaves={} # 记录已经计算过的原子最大值方向

class Calculator:
    def __init__(self,mol:Mol) -> None:
        self.mol=mol
        self.noNorms=[] # 已经计算过的原子法向量

    def pOrbital(self)->np.ndarray|None:
        """计算p轨道的方向向量"""
        pass
    
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
                wfn=wfnCaler.obtWfns(pos+delPos,[obt])[0]
                vals+=coef*wfn # 分子轨道=组合系数*原子轨道

            maxIdx=np.argmax(vals) # 最大值所在的索引
            maxVal=np.max(vals) # 最大值
            maxDir=delPos[maxIdx]
            return maxVal,maxDir
        from pywfn.spaceprop import wfnfunc
        wfnCaler=wfnfunc.Calculator(self.mol) # 波函数计算器
        step=0.1
        
        CM=self.mol.CM.copy()
        nmat=CM.shape[0]
        idxs=[] # 符合要求的轨道索引
        for i in range(nmat):
            if self.mol.obtAtms[i]!=atm:continue #不是指定的原子不算
            obtSym=self.mol.obtSyms[i]
            if re.match(sym,obtSym) is None:continue # 符号不对不算
            idxs.append(i)
        
            
        pos0=self.mol.atom(atm).coord.copy()
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
        
        direction=pos0-self.mol.atom(atm).coord
        length=np.linalg.norm(direction)
        if length<1e-8:
            return None
        else:
            return direction/length
        
    def sphAround(self)->np.ndarray|None: # 直接用dft的角度网格算了
        """计算周围一圈的球形范围"""
        from pywfn.data import lebedev
        lebedev.LD0074()
        

    def reactions(self,atm:int)->np.ndarray:
        """ reaction around
        计算原子可能的反应方向 以原子为中心的向量

        - 直线sp原子，垂直于键轴
        - 弯曲sp原子，键轴夹角垂面，且与键轴夹角小于90
        - 平面sp2，垂直于平面的两个方向
        - 弯曲sp2，垂直于平面，且与键轴之间夹角大于90°
        """
        def get_anyv(axis:np.ndarray): # 获取任意向量，但也不是真的任意
            for i in range(3):
                anyv=np.zeros(3)
                anyv[i]=1.0
                if vector_angle(anyv,axis)>1e-1:
                    return anyv
            raise ValueError("no anyv")
        atom=self.mol.atom(atm)
        nebs=atom.neighbors
        if len(nebs)==1:
            dirs=np.array([self.normal(atm)])
            return dirs
        if len(nebs)==2: # 两个原子的时候，代表卡宾C原子，在整个区域内分布
            ia,ib=nebs
            va=self.mol.atom(ia).coord-atom.coord
            vb=self.mol.atom(ib).coord-atom.coord
            va/=np.linalg.norm(va)
            vb/=np.linalg.norm(vb)
            vm=(va+vb)/2.0 # 中间向量

            angle=vector_angle(va,vb) # 两向量之间的夹角
            linear=abs(1-angle)<2e-2 # 是否为线性
            cent=atom.coord # 旋转中心，原子所在坐标
            # print('线性:',linear,abs(1-angle))
            
            if linear:
                axis=va-vb # 旋转轴
                anyv=get_anyv(axis) # 保证相同分子每次计算结果都一致
                cros=np.cross(axis,anyv) # 计算垂直于键轴和任意向量的向量
                cros/=np.linalg.norm(cros)
                angs=np.linspace(0,np.pi*2,35,endpoint=False)
                points=(cent+cros).reshape(-1,3)
                dirs=[points_rotate(points,cent,axis,ang)-cent for ang in angs]
                dirs=np.concatenate(dirs,axis=0)
            else:
                cros=np.cross(va,vb) # 垂直于va和vb的向量
                angs=np.linspace(0,np.pi,19,endpoint=True)
                cros/=np.linalg.norm(cros)
                points=(cent+cros).reshape(-1,3)
                dirs=[]
                for axis in (va,va-vb,-vb):
                    dirs+=[points_rotate(points,cent,axis,ang)-cent for ang in angs]
                dirs=np.concatenate(dirs,axis=0)
                if vector_angle(va+vb,dirs[28])<0.5:dirs*=-1
            return dirs
        if len(nebs)==3:
            pa,pb,pc=[self.mol.atom(n).coord for n in nebs]
            va,vb,vc=[p-atom.coord for p in (pa,pb,pc)]
            params=maths.get_plane_by_3points(pa,pb,pc)
            distan=maths.get_distance_to_plane(params,atom.coord)
            planer=distan<1e-1 # 是否为平面
            a,b,c,_=params
            normal=np.array([a,b,c])
            normal=normal/np.linalg.norm(normal)
            
            if planer:
                return np.array([normal,-normal])
            else:
                vmena=(va+vb+vc)/3.0
                angle=vector_angle(vmena,normal)
                if angle<0.5:normal*=-1
                vects=[normal]
                return np.array(vects)
        
        raise ValueError(f"原子{atm}未定义可能的反应方向!")
    
    def normal(self,atm:int)->np.ndarray|None:
        """计算原子的法向量

        Args:
            atm (int): 指定原子索引

        Returns:
            np.ndarray|None: 原子法向量，可能没有(None)
        """
        def deep_search(atm:int):
            atom=self.mol.atom(atm)
            searchd=[atm] # 已经搜索的原子
            searchs=atom.neighbors # 将要搜索的原子
            while len(searchs)>0:
                idx=searchs.pop(0)
                searchd.append(idx)
                atom=self.mol.atom(idx)
                norm=self.normal(idx)
                if norm is not None:
                    # print(idx,norm)
                    return idx,norm
                for each in atom.neighbors:
                    if each in searchd:continue
                    searchs.append(each)
            return None
        
        atom=self.mol.atom(atm)
        if 'normal' in atom._props.keys(): # 方便用户指定
            return atom._props['normal']
        nebs=atom.neighbors
        if len(nebs)==1: # 如果只连接一个原子
            searchs=deep_search(atm)
            assert searchs is not None,"未搜索到想要的原子"
            idx,normal=searchs
            return normal
        if len(nebs)==2:
            ia,ib=nebs
            va=self.mol.atom(ia).coord-atom.coord
            vb=self.mol.atom(ib).coord-atom.coord
            angle=vector_angle(va,vb)
            linear=abs(1-angle)<1e-1
            if linear:
                searchs=deep_search(atm)
                assert searchs is not None,"未搜索到想要的原子"
                idx,normal=searchs
                atom._props['normal']=normal
                return normal
            else:
                va/=np.linalg.norm(va)
                vb/=np.linalg.norm(vb)
                normal=np.cross(va,vb)
                normal/=np.linalg.norm(normal)
                if vector_angle(config.BASE_VECTOR,normal)>0.5:normal*=-1
                atom._props['normal']=normal
                return normal
        if len(nebs)==3:
            pa,pb,pc=[self.mol.atom(n).coord.copy() for n in nebs]
            vab=pb-pa
            vac=pc-pa
            vab/=np.linalg.norm(vab)
            vac/=np.linalg.norm(vac)
            normal=np.cross(vab,vac)
            normal/=np.linalg.norm(normal)
            if vector_angle(config.BASE_VECTOR,normal)>0.5:normal*=-1
            atom._props['normal']=normal
            return normal
        return None

    def coordSystem(self,atm:int,neb:int)->np.ndarray:
        """原子之上建立一组基坐标
        - x方向: atm -> beb
        - z方向: atm的法向量
        - y方向: y.z叉乘方向

        Args:
            atm (int): 要计算的原子索引
            neb (int): 相邻的原子索引

        Returns:
            np.ndarray: 坐标系,每一列代表一组基坐标
        """
        vx=self.mol.atom(neb).coord-self.mol.atom(atm).coord
        vx/=np.linalg.norm(vx)
        vz=self.normal(atm)
        vy=np.cross(vx,vz) # type: ignore
        return np.array([vx,vy,vz]).T
    
    def bases(self)->dict[int,np.ndarray]:
        """计算每个原子的局部基坐标[natm,3,3]，每一列代表一个基向量
        """
        atms=self.mol.heavyAtoms
        # base=np.zeros(shape=(natm,3,3))
        base={}
        for i,atm in enumerate(atms):
            atom=self.mol.atom(atm)
            norm=self.normal(atom.idx) # 原子法向量
            nebs=atom.neighbors
            
            cent=self.mol.coords[[e-1 for e in nebs],:].copy().mean(axis=0) # 邻居原子坐标的平均值
            vect=atom.coord-cent # 原子到邻居原子的向量
            vect/=np.linalg.norm(vect) # 单位化
            assert norm is not None,"原子法向量计算失败"
            if vector_angle(norm,vect)>0.5:norm*=-1
            anyv=np.random.rand(3) # 任意向量
            vz=norm/np.linalg.norm(norm) # 单位法向量
            vx=np.cross(vz,anyv)
            vx/=np.linalg.norm(vx)
            vy=np.cross(vz,vx)
            vy/=np.linalg.norm(vy)
            base[atm]=np.array([vx,vy,vz]).T # 原子局部坐标系的基坐标向量
        return base
    
    def hmoBases(self)->dict[int,np.ndarray]:
        """休克尔分子轨道法原子的方向"""
        from pywfn.orbtprop import popul
        NM=popul.Calculator(self.mol).piEleMatDecom('atom',False)
        eles=np.sum(NM,axis=0)
        piIdx=0 # 第一个pi分子轨道的索引
        for i,each in enumerate(eles):
            if abs(each-2)>0.5:continue
            piIdx=i
            break

        atms=self.mol.heavyAtoms
        base={}
        for i,atm in enumerate(atms):
            atom=self.mol.atom(atm)
            norm=self.normal(atom.idx) # 原子法向量
            if norm is None:continue
            u,l=atom.obtBorder
            coefs=self.mol.CM[u:l,piIdx]
            syms=self.mol.obtSyms[u:l]
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


def search_sp2Dir(v0,v1,v2,v3):
    """找到与va夹角为90°，且与vb,vc夹角相同的向量
    """
    def loss_fn(v):
        a1=vector_angle(v,v1)
        a2=vector_angle(v,v2)
        a3=vector_angle(v,v3)
        return (a1-0.5)**2+(a2-a3)**2
    def grad_fn(v):
        loss=loss_fn(v)
        grad=np.zeros(3)
        for i in range(3):
            dv=np.zeros(3)
            dv[i]=1e-3
            grad[i]=(loss_fn(v+dv)-loss)/1e-3
        return grad

    step=1.0
    loss0=np.inf
    for i in range(1000):
        grad=grad_fn(v0)
        vn=v0-step*grad
        vn/=np.linalg.norm(vn)
        loss=loss_fn(vn)
        # print(i,loss,step,grad)
        if loss<loss0:
            loss0=loss
            v0=vn
        else:
            if step<1e-6:
                break
            else:
                step/=10
    return v0