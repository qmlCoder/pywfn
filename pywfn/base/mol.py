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
from pywfn.maths import vector_angle,Gto
from pywfn import maths
from pywfn.base.atom import Atom,Atoms
from pywfn.base.bond import Bond,Bonds
from pywfn import reader
from pywfn import data
from pywfn.data.elements import elements
from pywfn import config
from pywfn import utils
from pywfn import reader

printer=utils.Printer()
import numpy as np
from functools import cached_property, lru_cache
from tqdm import tqdm
import multiprocessing as mp
from multiprocessing.pool import AsyncResult

class Mol:
    def __init__(self,reader:"reader.Reader"=None) -> None:
        self._atoms:Atoms=Atoms(self)
        self._bonds:Bonds=Bonds(self)
        self.reader:"reader.Reader"=reader
    
    @cached_property
    def basis(self)->data.Basis: # 设为属性，可以保证用不到的时候不被实例化
        return self.reader.get_basis()
    
    @cached_property
    def gto(self)->Gto:
        return Gto(self)


    @cached_property
    def charge(self):
        return self.reader.get_charge()
    
    @cached_property
    def spin(self):
        return self.reader.get_spin()

    @property
    def isOpenShell(self)->bool:
        """是否为开壳层"""
        CM=self.reader.get_CM()
        w,h=CM.shape
        return w!=h

    @cached_property
    def obtEcts(self):
        """每个分子轨道内的电子数量"""
        types=self.reader.get_obtTypes()
        oe=1 if self.isOpenShell else 2
        return [oe if o[-1]=='O' else 0 for o in types]

    @cached_property
    def obtEngs(self):
        return self.reader.get_obtEngs()

    @cached_property
    def obtStr(self)->list[str]:
        """返回轨道符号"""
        strs=[]
        for i,s in enumerate(self.reader.get_obtTypes()):
            if self.isOpenShell:
                if i<len(self.obtEcts)//2:
                    s=f'α {s}'
                else:
                    s=f'β {s}'
            idx=i
            if self.isOpenShell and i>=len(self.obtEcts)//2:
                idx=i-len(self.obtEcts)//2
            strs.append(f'{idx+1} {s}')
        return strs
    
    @cached_property
    def obtAtoms(self)->list[int]:
        return self.reader.get_obtAtoms()
    
    @cached_property
    def obtLayer(self)->list[str]:
        return self.reader.get_obtLayer()

    
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
                self._bonds.add(atom1,atom2)
        return self._bonds
    
    def atom(self,idx:int)->Atom:
        """根据原子编号获取一个原子，从1开始"""
        assert isinstance(idx,int),f"索引应该为整数"
        assert idx>0,"索引从1开始"
        return self.atoms[idx-1]
    
    @lru_cache
    def bond(self,idx1:int,idx2:int)->Bond:
        """第一次调用是生成键,第二次调用时直接返回,秒啊"""
        return self.bonds.get(idx1,idx2)

    @cached_property
    def coords(self)->np.ndarray:
        """返回原子坐标矩阵[n,3]"""
        return np.array([atom.coord for atom in self.atoms])

    @cached_property
    def O_obts(self)->list[int]:
        return [i for i,e in enumerate(self.obtEcts) if e!=0]
    
    @cached_property
    def V_obts(self)->list[int]:
        return [i for i,e in enumerate(self.obtEcts) if e==0]
    
    @cached_property
    def heavyAtoms(self)->list[Atom]:
        return [atom for atom in self.atoms if atom.symbol!='H']

    @cached_property
    def CM(self)->np.ndarray:
        """分子轨道系数矩阵"""
        return self.reader.get_CM()
    
    @property
    def SM(self):
        """重叠矩阵"""
        return self.reader.get_SM()
    
    @property
    def oE(self):
        """轨道电子数"""
        return 1 if self.isOpenShell else 2
    
    
    
    def projCM(self,atoms:list[int],obts:list[int],vects:list[np.ndarray],zero:bool=False):
        """
        获取投影后的系数矩阵
        atoms:需要投影的原子
        obts:需要投影的轨道
        vects:投影到的方向
        zero:其他系数是否置0
        atoms和vects的长度必须相同
        """
        assert isinstance(vects,list),"方向想两需要为列表"
        assert len(atoms)==len(vects),"原子的长度和方向的长度不同"
        
        if zero:
            CM_=np.zeros_like(self.CM,dtype=np.float32) #新的系数矩阵
        else:
            CM_=np.copy(self.CM)
        atoms:list[Atom]=(self.atom(a) for a in atoms)
        for a,(atom,vect) in enumerate(zip(atoms,vects)):

            a_1,a_2=atom.obtBorder
            layers=self.obtLayer[a_1:a_2]
            pIdx=[i for i,l in enumerate(layers) if 'P' in l]
            for o,obt in enumerate(obts):
                Co=np.zeros(len(layers)) #根据是否P轨道之外的保留还是0由不同的选择
                Cop=atom.get_pProj(vect,obt)
                Co[pIdx]=np.concatenate(Cop)
                CM_[a_1:a_2,obt]=Co
        return CM_

    def __repr__(self):
        return f'atom number: {len(self.atoms)}'

    def get_cloud(self,obt:int,pos:np.ndarray,atoms:list[int])->np.ndarray:
        """
        导出指定点的分子轨道数据
        obt: 要渲染的分子轨道
        pos: 要渲染的格点坐标[n,3]
        atoms: 要渲染的原子，索引从1开始
        """
        print(f'渲染轨道{obt=},{pos.shape},{atoms}')
        assert pos.shape[1]==3,'坐标的形状需为[n,3]'
        values=np.zeros(len(pos))
        pool=mp.Pool(mp.cpu_count())
        proces:list[AsyncResult]=[]
        # multi=True
        # if multi:
        #     for idx in atoms:
        #         atom=self.atom(idx)
        #         posi=pos-atom.coord
        #         proce=pool.apply_async(atom.get_cloud,(posi,obt,))
        #         proces.append(proce)
        #     for proce in tqdm(proces):
        #         values+=proce.get()
        #     return values
        # else:
        for idx in atoms:
            atom=self.atom(idx)
            posi=pos-atom.coord
            values+=atom.get_cloud(posi,obt)
        return values