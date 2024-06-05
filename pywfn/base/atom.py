from collections.abc import Iterator
import numpy as np
from functools import cached_property,lru_cache

from pywfn import config
from pywfn import base
from pywfn import maths

from pywfn.utils import printer

"""
一个原子的轨道组合系数就是一个矩阵，行数是基函数的数量，列数是分子轨道的数量
与轨道系数相关的太多了,尽量精简
OC是该原子的轨道系数
OCI为某一列的原子轨道系数,P:bool参数用来确定是否只要P轨道的系数
一个人原子的属性应该是标准的
原子轨道layer:List[str]
原子轨道系数OC:np.ndarray
读取器直接设置这些属性
"""
from pywfn.data.elements import elements
class Atom:
    def __init__(self,symbol:str,coord:np.ndarray,idx:int,mol:"base.Mol"): # 每个原子应该知道自己属于哪个分子
        self.symbol=symbol
        self.coord=coord
        self.coord.flags.writeable=False
        self.idx=idx
        self.mol:"base.Mol"=mol

        self._layersData={}
        self.layers:list[str]|None=None
        self._squareSum=None
        self._sContribution:dict={}
        self.OC:np.ndarray|None=None
        self._props:dict={}
    
    @cached_property
    def atomic(self)->int:
        return elements[self.symbol].idx
    
    @cached_property
    def radius(self)->float:
        return elements[self.symbol].radius

    @cached_property
    def obtBorder(self): # 获取元素上下界
        idxs=[i for i,idx in enumerate(self.mol.obtAtms) if idx==self.idx]
        return idxs[0],idxs[-1]+1 #因为最后一个元素不包含
    
    @property
    def obtCoeffs(self):
        """获取该原子对应的系数"""
        u,l=self.obtBorder
        return self.mol.CM[u:l,:]

    @cached_property
    def neighbors(self)->list[int]:
        """每个原子相邻的原子有哪些，根据分子的键来判断"""
        idxs=[]
        bonds=self.mol.bonds
        for bond in bonds:
            idx1,idx2=bond.ats
            if idx1==self.idx:
                idxs.append(idx2) # 添加不是该原子的键上的原子
            if idx2==self.idx:
                idxs.append(idx1)
        return [idx for idx in set(idxs)]
    

    @lru_cache
    def pLayersCs(self,obt:int):
        '''获取原子某一轨道的p层数据'''
        a,b=self.obtBorder
        layers=self.mol.obtSyms[a:b]
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
        a,b=self.obtBorder
        
        syms=self.mol.obtSyms[a:b]
        sidx=[i for i,l in enumerate(syms) if 'S' in l]
        return self.obtCoeffs[sidx,obt]
    
    
    # @lru_cache
    # def get_Normal(self,aroundIdx:int=None,main:bool=True)->None|np.ndarray: # 一个原子应该知道自己的法向量是多少
    #     """
    #     mian:代表是否为主动调用(防止递归)
    #     获取原子的法向量,垂直于键轴,如果不传入另一个原子,则垂直于三个原子确定的平面\n
    #     1.如果该原子连接两个原子
    #         1.1 如果两个原子不在一条直线上,有法向量，直接返回\n
    #         2.2 如果两个原子在一条之下上则没有法向量
    #     2.如果该原子没有法向量\n
    #         2.1 如果相邻的原子有法向量,返回邻原子的法向量\n
    #         2.2 如果相邻原子没有法向量,返回None
    #     """

    #     locNum=0.001
    #     neighbors=self.neighbors
    #     normal=None
    #     if len(neighbors)==3:
    #         if aroundIdx is None: #如果传入的没有相邻原子的序号
    #             #获取中心原子到三个相邻原子的法向量
    #             p1,p2,p3=[(each.coord-self.coord) for each in neighbors] 
    #         else:
    #             p1,p2,p3=[(each.coord-self.coord)*(1 if each.idx==aroundIdx else locNum) for each in neighbors] 
    #         normal=maths.get_normalVector(p1,p2,p3)
    #     elif len(neighbors)==2:
    #         p2=self.coord
    #         p1,p3=neighbors[0].coord,neighbors[1].coord
    #         angle=maths.vector_angle(p1-p2,p3-p2)
    #         if 0.2<=angle<=0.8:
    #             normal=maths.get_normalVector(p1,p2,p3,linear=True)
    #         else:
    #             pass
    #     elif len(neighbors)==1:
    #         a2=neighbors[0]
    #         normal = a2.get_Normal(self.idx)
    #     else:
    #         pass
    #     # elif main: #如果是在递归中调用本函数的话就不要再次递归了
    #     #     for each in neighbors:
    #     #         normal_=each.get_Normal(self.idx,main=False)
    #     #         if normal_ is not None:
    #     #             normal=normal_
    #     #             return normal
    #     #         return None #如果周围原子也没有法向量的话，返回None
    #     if normal is None:
    #         printer.warn(f'原子: {self.idx} 没有法向量')
    #     else:
    #         if maths.vector_angle(normal,config.BASE_VECTOR)>0.5:
    #             normal*=-1
    #     return normal


    # @lru_cache
    # def get_vertObt(self,idx1:int,idx2:int)->np.ndarray:
    #     """从HOMO轨道开始寻找垂直于键轴的轨道方向"""
    #     obts=self.mol.O_obts
    #     bondVector=self.mol.atom(idx1).coord-self.mol.atom(idx2).coord
    #     for obt in obts[::-1]:
    #         obtWay=self.get_obtWay(obt)
    #         if obtWay is None:continue
    #         printer.info(f'以{obt}号轨道方向作为法向量')
    #         if maths.vector_angle(bondVector,obtWay,trans=True)>0.4:
    #             if maths.vector_angle(obtWay,config.BASE_VECTOR)>0.5:obtWay*=-1
    #             return obtWay
    
    # def get_projWay(self):
    #     """
    #     获得系数投影方向
    #     有法向量用法向量,没有法向量用垂直于轨道和键轴的方向
    #     没有法向量的时候,相当于原子在线性区域,两个键平行,但还是用平均值比较合理
    #     """
    #     norm=self.get_Normal()
    #     if norm is not None:
    #         return norm
    #     else:
    #         a1,a2=self.neighbors #相邻两个原子
    #         return self.get_vertObt(a1,a2) #返回垂直于键轴的轨道方向
    
    # @lru_cache
    # def get_obtWay(self,obt:int)->np.ndarray:
    #     """获得原子某一轨道的方向"""
    #     maxPos,maxValue=maths.get_extraValue(self,obt)
    #     way=maxPos # 如果两者相同说明完全没有电子云、这显然不对啊
    #     if np.linalg.norm(way)==0:
    #         printer.warn(f'轨道{obt}无方向')
    #         return None
    #     return way/np.linalg.norm(way)

    # def __repr__(self) -> str:
    #     x,y,z=self.coord
    #     return f'{self.idx} {self.symbol} ({x:.4f},{y:.4f},{z:.4f})'

    # @lru_cache
    # def get_sCont(self,orbital:int):
    #     """获取某个原子轨道的贡献"""
    #     s=self.OC[0,orbital]
    #     contribution=s**2/self.mol.As[orbital]
    #     return contribution

    def get_direct(self):
        """
        获取反应方向
        对于两个键的，为平面内夹角相同方向
        
        """
    
    def __repr__(self) -> str:
        return f'Atom:({self.symbol},{self.idx})'
        


class Atoms:
    def __init__(self,mol:"base.Mol") -> None:
        self.mol=mol
        self.atoms:list[Atom]=[]
    
    def __repr__(self) -> str:
        atoms=[atom.symbol for atom in self.atoms]
        return f'Atoms:{len(self.atoms)}\n'+' '.join(atoms)
    
    def add(self,symbol:str,coord:np.ndarray):
        atom=Atom(symbol,coord,len(self.atoms)+1,self.mol)
        self.atoms.append(atom)
    
    @property
    def atomics(self)->list[int]:
        return [a.atomic for a in self.atoms]
    
    @property
    def symbols(self)->list[str]:
        return [a.symbol for a in self.atoms]
    
    
    @property
    def indexs(self)->list[int]:
        return [a.idx for a in self.atoms]
    
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