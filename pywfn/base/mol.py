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
from tqdm import tqdm
import multiprocessing as mp
from multiprocessing.pool import AsyncResult
import threading
from typing import Callable

class Props(dict):
    def __init__(self) -> None:
        super().__init__({})

    def get(self,key:str,fun:Callable):
        if key not in self.keys():
            self[key]=fun()
        return self[key]

    def set(self,key:str,value):
        self[key]=value

class Mol:
    """基础的分子对象"""
    def __init__(self,reader:"reader.Reader"=None) -> None:
        """
        分子实例化
        `reader:Reader`，分子读取器
        """
        self._atoms:Atoms=Atoms(self)
        self._bonds:Bonds=Bonds(self)
        self.reader:"reader.Reader"=reader
        self.props=Props()
        self.bohr:bool=False # 是否使用波尔坐标
    
    @cached_property
    def basis(self)->data.Basis: # 设为属性，可以保证用不到的时候不被实例化
        return self.reader.get_basis()
    
    @cached_property
    def gto(self)->"maths.Gto":
        return maths.Gto(self)

    @cached_property
    def charge(self)->int:
        return self.reader.get_charge()
    
    @cached_property
    def spin(self)->int:
        return self.reader.get_spin()

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
                if i<len(nobt)//2:
                    s=f'α {s}'
                else:
                    s=f'β {s}'
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

    @property
    def atoms(self)->Atoms:
        """获取所有原子"""
        if self._atoms:return self._atoms
        symbols=self.reader.get_symbols()
        coords=self.reader.get_coords()
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
                r_=elements[atom1.symbol].radius+elements[atom2.symbol].radius
                if r>r_*1.1:continue
                self._bonds.add(atom1.idx,atom2.idx)
        return self._bonds
    
    def atom(self,idx:int)->Atom:
        """根据原子编号获取一个原子，从1开始"""
        assert isinstance(idx,int),f"索引应该为整数"
        assert idx>0,"索引从1开始"
        assert idx<=len(self.atoms),"原子索引超过原子数量"
        return self.atoms[idx-1]
    
    @lru_cache
    def bond(self,idx1:int,idx2:int)->Bond:
        """第一次调用是生成键,第二次调用时直接返回,秒啊"""
        return self.bonds.get(idx1,idx2)

    @cached_property
    def coords(self)->np.ndarray:
        """返回原子坐标矩阵[n,3]"""
        return np.array([atom.coord for atom in self.atoms])

    @property
    def O_obts(self)->list[int]:
        return [i for i,occ in enumerate(self.obtOccs) if occ]
    
    @cached_property
    def V_obts(self)->list[int]:
        return [i for i,occ in enumerate(self.obtOccs) if not occ]
    
    @cached_property
    def heavyAtoms(self)->list[Atom]:
        return [atom for atom in self.atoms if atom.symbol!='H']

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
        return maths.CM2PM(self.CM,self.O_obts,self.oE)
    
    @property
    def oE(self):
        """轨道电子数"""
        return 1 if self.open else 2
    
    def projCM(self,obts:list[int],atms:list[int],dirs:list[np.ndarray]
               ,akeep:bool,lkeep:bool,ins:bool)->np.ndarray:
        """
        获取投影后的系数矩阵
        atms:需要投影的原子
        obts:需要投影的轨道
        dirs:投影到的方向 atms和dirs的长度必须相同
        akeep:其它原子系数是否保留 keep other atoms
        lkeep:其它价层系数是否保留 keep other layer
        ins:是否包含价层s轨道
        """
        assert isinstance(dirs,list),"方向想两需要为列表"
        assert len(atms)==len(dirs),"原子和方向数量不同"
        if akeep:
            CMp=np.copy(self.CM)
        else:
            CMp=np.zeros_like(self.CM,dtype=np.float32) #新的系数矩阵
            
        for a,(atom,vect) in enumerate(zip(atms,dirs)):
            atom=self.atom(atom)
            nebNum=len(atom.neighbors)
            u,l=atom.obtBorder
            syms=self.obtSyms[u:l]
            
            pIdx=[i for i,s in enumerate(syms) if 'P' in s]
            if len(pIdx)==0:continue # 没有p轨道则跳过
            for o,obt in enumerate(obts):
                if lkeep:
                    Co=self.CM.copy()[u:l,obt]
                else:
                    Co=np.zeros(len(syms)) #根据是否P轨道之外的保留还是0由不同的选择
                if ins:
                    sIdx=[i for i,s in enumerate(syms) if 'S' in s][1:]
                    Cos=atom.obtCoeffs.copy()[sIdx,obt]
                    Co[sIdx]=Cos*(3-nebNum)/3
                Cop=atom.get_pProj(vect,obt)
                Co[pIdx]=np.concatenate(Cop)
                CMp[u:l,obt]=Co.copy()
        return CMp

    def __repr__(self):
        return f'atom number: {len(self.atoms)}'

    def get_wfnv(self,obt:int,pos:np.ndarray,atoms:list[int])->np.ndarray:
        """
        导出指定点的分子轨道数据
        obt: 要渲染的分子轨道
        pos: 要渲染的格点坐标[n,3]
        atoms: 要渲染的原子，索引从1开始
        """
        print(f'渲染轨道{obt=},{pos.shape},{atoms}')
        assert pos.shape[1]==3,'坐标的形状需为[n,3]'
        values=np.zeros(len(pos))
        # pool=mp.Pool(mp.cpu_count())
        # proces:list[AsyncResult]=[]
        # multi=True
        # if multi:
        #     for idx in atoms:
        #         atom=self.atom(idx)
        #         posi=pos-atom.coord
        #         proce=pool.apply_async(atom.get_wfnv,(posi,obt,))
        #         proces.append(proce)
        #     for proce in tqdm(proces):
        #         values+=proce.get()
        #     return values
        # else:
        for idx in atoms:
            atom=self.atom(idx)
            posi=pos-atom.coord
            values+=atom.get_wfnv(posi,obt)
        return values
    
    def get_dens(self,atoms:list[int],obts:list[int],coords:np.ndarray)->np.ndarray:
        """
        计算分子的电子密度
        atoms:要计算的原子，从1开始
        coords:空间笛卡尔坐标
        """
        molDens=np.zeros(shape=(len(coords))) # 每一个坐标的电子密度值
        for a in atoms:
            atom=self.atom(a)
            coord=coords-atom.coord #以原子为中心的坐标
            atmDens=atom.get_dens(obts,coord)
            molDens+=atmDens
        return molDens