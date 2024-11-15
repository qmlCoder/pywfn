"""
基础的分子对象，其属性应该是标准的，已知的
一个分子对象应该有哪些属性？基本属性(必须属性),计算属性(需要计算才能的到的属性)
结构相关
- 原子
- 键
- 法向量
- 重原子

分子轨道相关
轨道数量
轨道类型
是否分为α和β
占据轨道序数
非占据轨道序数
重叠矩阵

读取器
"""

from typing import Any
from pywfn import maths
from pywfn.base.atom import Atom,Atoms
from pywfn.base.bond import Bond,Bonds
from pywfn import reader
from pywfn import data
from pywfn.data.elements import elements
from pywfn import config
from pywfn import utils
from pywfn import reader
from pywfn.utils import printer

import numpy as np
from functools import cached_property, lru_cache
import multiprocessing as mp
from multiprocessing.pool import AsyncResult
import threading
from typing import Callable
import collections
import re
import copy

class Props(dict):
    def __init__(self) -> None:
        super().__init__({})

    def get(self,key:str,fun:Callable):
        if key not in self.keys():
            self.set(key,fun())
        return self[key]

    def set(self,key:str,value):
        self[key]=value

class Mol:
    """基础的分子对象"""
    def __init__(self,reader:"reader.Reader") -> None:
        """
        分子实例化
        `reader:Reader`，分子读取器
        """
        self._atoms:Atoms=Atoms(self)
        self._bonds:Bonds=Bonds(self)
        self.reader:"reader.Reader"=reader
        self.props:Props=Props()
        self.bohr:bool=False # 是否使用波尔坐标
    
    @cached_property
    def basis(self)->data.Basis: # 设为属性，可以保证用不到的时候不被实例化
        return self.reader.get_basis()
    
    @cached_property
    def gto(self)->"maths.Gto":
        return maths.Gto(self)

    @property
    def charge(self)->int:
        charge=self.props.get('charge',self.reader.get_charge)
        return charge
    
    @property
    def spin(self)->int:
        spin=self.props.get('spin',self.reader.get_spin)
        return spin

    @property
    def open(self)->bool:
        """是否为开壳层"""
        w,h=self.CM.shape
        return w!=h
    
    @property
    def energy(self)->float:
        """获取分子能量"""
        return self.props.get('energy',self.reader.get_energy)

    @property
    def obtOccs(self)->list[bool]:
        """获取每个分子轨道是否占据"""
        occs=self.props.get('obtOccs',self.reader.get_obtOccs)
        return occs

    @cached_property
    def obtEngs(self)->list[float]:
        return self.reader.get_obtEngs()

    @cached_property
    def obtStrs(self)->list[str]:
        """返回轨道符号"""
        strs=[]
        occs=self.obtOccs
        nobt=len(occs)
        for i,s in enumerate(occs):
            if self.open:
                if i<nobt//2:
                    s=f'A {s}'
                else:
                    s=f'B {s}'
            idx=i
            if self.open and i>=nobt//2:
                idx=i-nobt//2
            strs.append(f'{idx+1} {s}')
        return strs
    
    @cached_property
    def obtAtms(self)->list[int]:
        return self.reader.get_obtAtms()
    
    @cached_property
    def obtShls(self)->list[int]:
        return self.reader.get_obtShls()

    @cached_property
    def obtSyms(self)->list[str]:
        return self.reader.get_obtSyms()
    
    @cached_property
    def obtAngs(self)->list[int]:
        syms=self.obtSyms
        angs=[self.basis.sym2ang(sym) for sym in syms]
        return angs

    @property
    def atoms(self)->Atoms:
        """获取所有原子"""
        if self._atoms:return self._atoms
        symbols=self.props.get('symbols',self.reader.get_symbols)
        coords=self.props.get('coords',self.reader.get_coords)
        for s,c in zip(symbols,coords):
            self._atoms.add(symbol=s,coord=c)
        return self._atoms
    
    @property
    def bonds(self)->Bonds:
        if self._bonds:return self._bonds
        for atom1 in self.atoms:
            for atom2 in self.atoms:
                if atom1.idx>=atom2.idx:continue
                r=np.linalg.norm(atom2.coord-atom1.coord)
                r1=elements[atom1.symbol].radius
                r2=elements[atom2.symbol].radius
                # if r>config.BOND_LIMIT:continue
                if r>(r1+r2)*1.1:continue
                self._bonds.add(atom1.idx,atom2.idx)
        return self._bonds
    
    def atom(self,idx:int)->Atom:
        """根据原子编号获取一个原子，从1开始"""
        assert isinstance(idx,int),f"索引应该为整数"
        assert idx>0,"索引从1开始"
        assert idx<=len(self.atoms),"原子索引超过原子数量"
        return self.atoms[idx-1]
    
    @lru_cache
    def bond(self,idx1:int,idx2:int)->Bond|None:
        """第一次调用是生成键,第二次调用时直接返回,秒啊"""
        return self.bonds.get(idx1,idx2)

    @property
    def coords(self)->np.ndarray:
        """返回原子坐标矩阵[n,3]"""
        coords=[atom.coord for atom in self.atoms]
        coords=np.array(coords)
        coords.setflags(write=False)
        return coords

    @property
    def O_obts(self)->list[int]:
        return [i for i,occ in enumerate(self.obtOccs) if occ]
    
    @cached_property
    def V_obts(self)->list[int]:
        return [i for i,occ in enumerate(self.obtOccs) if not occ]
    
    @cached_property
    def heavyAtoms(self)->list[int]:
        return [atom.idx for atom in self.atoms if atom.symbol!='H']

    @property
    def CM(self)->np.ndarray:
        """分子轨道系数矩阵"""
        return self.props.get('CM',self.reader.get_CM)
    
    @property
    def SM(self):
        """重叠矩阵"""
        return self.reader.get_SM()
    
    @property
    def PM(self):
        """密度矩阵"""
        return maths.CM2PM(self.CM.copy(),self.O_obts,self.oE)
    
    @property
    def oE(self):
        """轨道电子数"""
        return 1 if self.open else 2
    
    @property
    def formula(self)->str:
        """获取分子式"""
        symbols=self.atoms.symbols
        # 统计列表中每个元素出现的次数
        counts = collections.Counter(symbols)
        names=[f'{k}{v}' for k,v in counts.items()]
        return ''.join(names)
    
    @property
    def molBorder(self):
        p0=self.coords.min(axis=0)
        p1=self.coords.max(axis=0)
        return p0,p1
    

    def params(self,atms:tuple[int])->float:
        """
        获取键长、键角、二面角
        """
        if len(atms)==2:
            a1,a2=atms
            bond=self.atom(a1).coord-self.atom(a2).coord
            return float(np.linalg.norm(bond))
        elif len(atms)==3:
            a1,a2,a3=atms
            v21=self.atom(a1).coord-self.atom(a2).coord
            v23=self.atom(a3).coord-self.atom(a2).coord
            angle=maths.vector_angle(v21,v23)
            return angle
        elif len(atms)==4:
            a1,a2,a3,a4=atms
            v21=self.atom(a1).coord-self.atom(a2).coord
            v23=self.atom(a3).coord-self.atom(a2).coord
            v32=self.atom(a2).coord-self.atom(a3).coord
            v34=self.atom(a3).coord-self.atom(a4).coord
            n2=np.cross(v21,v23)
            n3=np.cross(v34,v32)
            n2/=np.linalg.norm(n2)
            n3/=np.linalg.norm(n3)
            angle=maths.vector_angle(n2,n3)
            nm=np.dot(v21,n3) # 为了实现二面角的正负而引入
            nm/=np.linalg.norm(nm)
            return float(angle*np.pi*nm)
        else:
            raise ValueError("参数数量错误")
    
    @property
    def eleNum(self)->tuple[int,int]:
        """获取alpha,beta电子数"""
        obts=self.O_obts
        if self.open:
            nbas=self.CM.shape[0]
            na=sum([1 for o in obts if o<nbas])
            nb=len(obts)-na
            return na,nb
        else:
            return len(obts),len(obts)
    
    @cached_property
    def DM(self):
        """计算键长矩阵"""
        natm=self.atoms.num
        DM=np.zeros(shape=(natm,natm))
        for i in range(natm):
            for j in range(i,natm):
                dis=self.params((i+1,j+1))
                DM[i,j]=dis
                DM[j,i]=dis
        return DM

    def __repr__(self):
        atmDict={}
        for atom in self.atoms:
            if atom.symbol not in atmDict:
                atmDict[atom.symbol]=1
            else:
                atmDict[atom.symbol]+=1
        atmList=[(k,v) for k,v in atmDict.items()]
        # print('排序前',atmList)
        atmList.sort(key=lambda x:x[1],reverse=True)
        # print('排序后',atmList)
        name=''.join([f'{k}{v}' for k,v in atmList])
        natm=self.atoms.num
        nele=sum(self.eleNum)
        return f'<{name},{natm},{nele},{self.charge},{self.spin}>'
    
    def clone(self):
        return copy.deepcopy(self)