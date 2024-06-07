"""
计算原子的各种方向
- 原子法向量
- p轨道方向
- 原子周围一圈方向
"""
from pywfn.base import Mol
from pywfn import maths
from pywfn.maths import vector_angle,points_rotate
from pywfn import config

import numpy as np

class Calculator:
    def __init__(self,mol:Mol) -> None:
        self.mol=mol

    def pObt(self)->np.ndarray|None:
        """计算p轨道的方向向量"""
        pass

    def maxWfn(self)->np.ndarray|None:
        """计算波函数最大值方向"""

    def sphAro(self)->np.ndarray|None:
        """计算周围一圈的球形范围"""
        pass

    def reaction(self,atm:int)->np.ndarray:
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
        if len(nebs)==2:
            ia,ib=nebs
            va=self.mol.atom(ia).coord-atom.coord
            vb=self.mol.atom(ib).coord-atom.coord
            va/=np.linalg.norm(va)
            vb/=np.linalg.norm(vb)
            vm=(va+vb)/2.0 # 中间向量

            angle=vector_angle(va,vb) # 两向量之间的夹角
            linear=abs(1-angle)<1e-1 # 是否为线性
            cent=atom.coord
            axis=va-vb # 要保证旋转
            if linear:
                anyv=get_anyv(axis) # 保证相同分子每次计算结果都一致
                cros=np.cross(axis,anyv) # 计算垂直于键轴和任意向量的向量
                angs=np.linspace(0,np.pi*2,35,endpoint=False)
            else:
                cros=np.cross(va,vb) # 垂直于va和vb的向量
                angs=np.linspace(0,np.pi,19,endpoint=True)
            cros/=np.linalg.norm(cros)
            points=(cent+cros).reshape(-1,3)
            dirs=[points_rotate(points,cent,axis,ang)-cent for ang in angs]
            dirs=np.concatenate(dirs,axis=0)
            if vector_angle(vm,dirs[9])>0.5:dirs*=-1
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
            vmena=(va+vb+vc)/3.0
            angle=vector_angle(vmena,normal)
            if angle<0.5:normal*=-1
            if planer:
                return np.array([normal,-normal])
            else:
                vects=[
                    normal,
                    search_sp2Dir(normal.copy(),va,vb,vc),
                    search_sp2Dir(normal.copy(),vb,va,vc),
                    search_sp2Dir(normal.copy(),vc,va,vb),
                ]
                return np.array(vects)
        
        raise ValueError("找不到可能的反应方向")
    
    def normal(self,atm:int)->np.ndarray|None:
        """计算原子的法向量"""
        
        atom=self.mol.atom(atm)
        if 'normal' in atom._props.keys(): # 方便用户指定
            return atom._props['normal']
        nebs=atom.neighbors
        if len(nebs)==1:
            neb=nebs[0]
            nebb=self.mol.atom(neb).neighbors
            if atom.symbol=='H': # H肯定没有法向量吧
                return None
            elif len(nebb) in [2,3]:
                return self.normal(neb)
            else:
                return None
        if len(nebs)==2:
            ia,ib=nebs
            va=self.mol.atom(ia).coord-atom.coord
            vb=self.mol.atom(ib).coord-atom.coord
            angle=vector_angle(va,vb)
            linear=abs(1-angle)<1e-1
            if linear:
                return None
            else:
                va/=np.linalg.norm(va)
                vb/=np.linalg.norm(vb)
                normal=np.cross(va,vb)
                normal/=np.linalg.norm(normal)
                if vector_angle(config.BASE_VECTOR,normal)>0.5:normal*=-1
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
            return normal
        return None

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