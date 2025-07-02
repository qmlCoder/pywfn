
from collections.abc import Iterator
import numpy as np
from functools import cached_property,lru_cache
from pywfn import config
from pywfn import base
from pywfn import maths
from pywfn.maths import vector_angle
from pywfn.utils import printer
from pywfn.data.elements import elements
class Atom:
    """
    原子对象
    """
    def __init__(self,sym:str,xyz:np.ndarray,idx:int,geome:"base.Geome"): # 每个原子应该知道自己属于哪个分子
        self.sym=sym
        self.xyz=xyz
        self.xyz.flags.writeable=False
        self.idx=idx
        self.geome=geome

        self._layersData={}
        self.layers:list[str]|None=None
        self._squareSum=None
        self._sContribution:dict={}
        self._props:dict={}
    
    @property
    def mol(self):
        assert self.geome.mol is not None,"未绑定分子"
        return self.geome.mol
    
    @property
    def symbol(self)->str:
        return self.sym
    
    @property
    def coord(self)->np.ndarray:
        return self.xyz

    @cached_property
    def atomic(self)->int:
        return elements[self.sym].idx
    
    @cached_property
    def radius(self)->float:
        return elements[self.sym].radius

    @cached_property
    def obtBorder(self): # 获取元素上下界
        idxs=[i for i,idx in enumerate(self.mol.atoAtms) if idx==self.idx]
        return idxs[0],idxs[-1]+1 #因为最后一个元素不包含
    
    @cached_property
    def obtCoeffs(self):
        """获取该原子对应的系数"""
        u,l=self.obtBorder
        return self.mol.CM[u:l,:]
    
    @cached_property
    def OC(self):
        return self.obtCoeffs

    @cached_property
    def neighbors(self)->tuple[int,...]:
        """每个原子相邻的原子有哪些，根据分子的键来判断"""
        idxs=[]
        bonds=self.mol.bonds
        for bond in bonds:
            idx1,idx2=bond.ats
            if idx1==self.idx:
                idxs.append(idx2) # 添加不是该原子的键上的原子
            if idx2==self.idx:
                idxs.append(idx1)
        return tuple(idxs)

    @lru_cache
    def pLayersCs(self,obt:int):
        '''获取原子某一轨道的p层数据'''
        u,l=self.obtBorder
        layers=self.mol.atoSyms[u:l]
        pIndex=[i for i,l in enumerate(layers) if 'P' in l]
        return self.obtCoeffs[pIndex,obt]
    
    def get_pProj(self,direct:np.ndarray,obt:int)->list[np.ndarray]:
        """计算原子p系数在某个方向上的投影,返回n个三维向量"""
        assert isinstance(direct,np.ndarray),"方向需要为np.ndarray"
        assert direct.shape==(3,),'方向向量长度应该为3'
        length=np.linalg.norm(direct)
        assert np.abs(length-1)<1e-4,'方向范数应该为1'
        assert length!=0,"方向向量长度不能为0"
        cs=self.pLayersCs(obt) # p轨道的系数
        assert len(cs)%3==0,"p轨道系数的长度应为3n"
        ps=[np.array(cs[i:i+3]) for i in range(0,len(cs),3)] #每一项都是长度为3的数组
        lens=np.dot(ps,direct)
        ps_=[l*direct for l in lens] # 轨道向量在法向量方向上的投影
        return ps_
    
    def get_sProj(self,direct:np.ndarray,obt:int):
        """计算原子p系数在某个方向上的投影,返回n个三维向量"""
        assert isinstance(direct,np.ndarray),"方向需要为np.ndarray"
        assert np.linalg.norm(direct)!=0,"方向向量长度不能为0"
        direct/=np.linalg.norm(direct) # 投影向量归一化
        u,l=self.obtBorder
        syms=self.mol.atoSyms[u:l]
        sidx=[i for i,l in enumerate(syms) if 'S' in l]
        return self.obtCoeffs[sidx,obt]
    
    @cached_property
    def is_linear(self):
        """判断原子是否是线性的"""
        if len(self.neighbors)!=2:return False
        a1,a2=list(self.neighbors)
        v1=self.mol.atom(a1).coord-self.coord
        v2=self.mol.atom(a2).coord-self.coord
        angle=vector_angle(v1,v2)
        return abs(1-angle)<1e-1

    
    def __repr__(self) -> str:
        return f'Atom:({self.symbol},{self.idx})'

class Atoms:
    def __init__(self,geome:"base.Geome") -> None:
        self.geome=geome
        self.atoms:list[Atom]=[]
    
    def __repr__(self) -> str:
        atoms=[atom.symbol for atom in self.atoms]
        return f'Atoms:{len(self.atoms)}\n'+' '.join(atoms)
    
    def add(self,symbol:str,coord:np.ndarray):
        assert coord.shape==(3,), '坐标必须是三维'
        atom=Atom(symbol,coord,len(self.atoms)+1,self.geome)
        self.atoms.append(atom)
    
    @property
    def atomics(self)->list[int]:
        return [a.atomic for a in self.atoms]
    
    @property
    def symbols(self)->list[str]:
        return [a.symbol for a in self.atoms]
    
    @property
    def syms(self)->list[str]:
        return [a.symbol for a in self.atoms]
    
    @property
    def xyzs(self)->np.ndarray:
        return np.array([a.coord for a in self.atoms])
    
    @property
    def indexs(self)->list[int]:
        return [a.idx for a in self.atoms]
    
    @property
    def atms(self)->list[int]:
        return [a.idx for a in self.atoms]
    
    @property
    def natm(self)->int:
        return len(self.atms)
    
    @cached_property
    def LM(self):
        """原子之间的键长矩阵"""
        nmat=self.num
        LM=np.zeros((nmat,nmat),dtype=np.float32)
        for i,atomi in enumerate(self.atoms):
            for j,atomj in enumerate(self.atoms):
                LM[i,j]=np.linalg.norm(atomj.coord-atomi.coord)
        return LM
    
    @property
    def radius(self):
        """所有原子的半径"""
        radius=[elements[atom.symbol].radius for  atom in self.atoms]
        return radius

    @property
    def atmuls(self)->list[tuple[int,int]]:
        """所有原子在基函数矩阵中的上下限"""
        atmuls=[atom.obtBorder for atom in self.atoms]
        return atmuls
    
    @property
    def uls(self)->list[tuple[int,int]]:
        """所有原子在基函数矩阵中的上下限"""
        atmuls=[atom.obtBorder for atom in self.atoms]
        return atmuls
    
    @property
    def num(self)->int:
        return len(self.atoms)

    def __getitem__(self,item)->Atom:
        return self.atoms[item]
    
    def __iter__(self)->Iterator[Atom]:
        for atom in self.atoms:
            yield atom
        
    def __bool__(self)->bool:
        return len(self.atoms)!=0
    
    def __len__(self)->int:
        return len(self.atoms)