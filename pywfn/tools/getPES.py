"""
本脚本用以生成势能面
"""
from pathlib import Path

from pywfn.base import Mole
from pywfn.reader import get_reader,LogReader
from pywfn.data import temps
from collections import namedtuple
from dataclasses import dataclass
from typing import Union

@dataclass
class Mole:
    energy:float # 分子能量
    name:str='' # 分子名

@dataclass
class Moles:
    mols:list[Mole]
    text:str='mol'
    n0:Union[None,"Node"]=None # 开始的节点
    n1:Union[None,"Node"]=None # 结束的节点
    posY:float=0
    @property
    def energy(self):
        return sum([mol.energy for mol in self.mols])
    

@dataclass
class Block:
    blocks:"Blocks"
    idx:int
    Rms:Moles # 反应物
    Tms:Moles # 过渡态
    Pms:Moles # 产物
    bi0:Union[int,"Block"]=None
    bi1:Union[int,"Block"]=None

    @property
    def b0(self)->"Block":
        if self.bi0 is None:
            return None
        else:
            return self.blocks[self.bi0]

    @property
    def b1(self)->"Block":
        if self.bi1 is None:
            return None
        else:
            return self.blocks[self.bi1]

    @property
    def engs(self)->list[float]:
        e0=self.Rms.energy
        es=[self.Rms.energy,self.Tms.energy,self.Pms.energy]
        es=[e-e0 for e in es] # 反应物为0的能量
        if self.b0:
            es=[e+self.b0.Pms.energy for e in es] # 接着上一块的能量
        return es
    
    @property
    def pos(self):
        if self.b0:
            return self.b0.pos+Sizes.block
        else:
            return 0
    
    def __repr__(self):
        return f'{self.bi0}<==({self.idx})==>{self.bi1}'

@classmethod
class Blocks:
    blocks:list[Block]=[]

    def __getitem__(self,idx)->Block:
        for block in self.blocks:
            if idx==block.idx:
                return block
    
    def append(self,block:Block):
        self.blocks.append(block)

@dataclass
class Node: # chemdraw中的一个节点
    id:int
    x:float
    y:float

@dataclass
class Bond:
    id:int
    B:int
    E:int
    Display:str

@dataclass
class Text:
    id:int
    t:str
    x:float
    y:float

class Sizes:
    block:float=90 # 一个块的大小
    moles:float=30 # 一个分子组的大小
    startx:float=50
    starty:float=100


class Tool:
    def __init__(self) -> None:
        self.blocks:list[Block]=[]
        # self.route:list[int]
        self.temp=temps.pes
        self.nodes:list[Node]=[]
        self.bonds:list[Bond]=[]
        self.texts:list[Text]=[]
        self.NBS=''
        self.TXS=''
        self.idx=4
    
    def getIdx(self):
        self.idx+=1
        return self.idx
    
    def getPos(self,x:float,y:float):
        return x+Sizes.startx+Sizes.block/2+10,400-(y+Sizes.starty)

    def create(self):
        for b,block in enumerate(self.blocks):
            self.add_block(block)
        self.build_bonds()
        self.write()

    def add_block(self,block:Block):
        """
        添加一个Block块
        bid:第多少个块
        """
        Rms=block.Rms
        Tms=block.Tms
        Pms=block.Pms
        eng0=0 if block.b0 is None else block.b0.Pms.posY
        eng1=eng0+(Tms.energy-Rms.energy)
        eng2=eng0+(Pms.energy-Rms.energy)
        Rms.posY,Tms.posY,Pms.posY=eng0,eng1,eng2
        
        
        x1=block.pos
        x2=x1+Sizes.block/2
        
        if block.b0 is None:
            x0=block.pos-Sizes.block/2
            Rms.n0,Rms.n1=self.add_moles(eng0,x0,Rms.text)
            
        Tms.n0,Tms.n1=self.add_moles(eng1,x1,Tms.text)
        Pms.n0,Pms.n1=self.add_moles(eng2,x2,Pms.text)
        print(block.idx,x1)
        

    def add_moles(self,eng:float,x:float,text:str): # 一个moles是一个横线
        x0=x
        x1=x+Sizes.moles
        n0=self.add_node(x0,eng)
        n1=self.add_node(x1,eng)
        self.add_text(text,(x0+x1)/2-len(text)*4,eng+2)
        return n0,n1
            
    def add_node(self,x,y)->Node:
        """添加一个节点"""
        x,y=self.getPos(x,y)
        node=Node(id=self.getIdx(),x=x,y=y)
        self.nodes.append(node)
        return node

    def build_bonds(self):
        """生成键"""
        for block in self.blocks:
            if block.b0 is None:
                self.add_bond(block.Rms.n0,block.Rms.n1,'Bold')
                self.add_bond(block.Rms.n1,block.Tms.n0,'Dash')
            if block.b0:
                self.add_bond(block.b0.Pms.n1,block.Tms.n0,'Dash')
            self.add_bond(block.Tms.n0,block.Tms.n1,'Bold')
            self.add_bond(block.Tms.n1,block.Pms.n0,'Dash')
            self.add_bond(block.Pms.n0,block.Pms.n1,'Bold')
            
    def add_bond(self,nb:Node,ne:Node,display:str):
        """添加键"""
        bond=Bond(id=self.getIdx(),B=nb.id,E=ne.id,Display=display)
        self.bonds.append(bond)
        return self.idx
    
    def add_text(self,text:str,x:float,y:float):
        """添加文本"""
        x,y=self.getPos(x,y)
        self.texts.append(Text(id=self.getIdx(),t=text,x=x,y=y))
        
    
    def write(self):
        for node in self.nodes:
            self.NBS+=f'<n id="{node.id}" p="{node.x} {node.y}" Z="1" AS="N"/>\n'
        for bond in self.bonds:
            self.NBS+=f'<b id="{bond.id}" B="{bond.B}" E="{bond.E}" Z="2" Display="{bond.Display}" BS="N"/>\n'
        for text in self.texts:
            self.TXS+=f'<t id="{text.id}" p="{text.x} {text.y}" Z="3" LineHeight="auto"><s font="5" size="10" color="0">{text.t}</s></t>\n'
        self.temp=self.temp.replace('[NBS]',self.NBS)
        self.temp=self.temp.replace('[TXS]',self.TXS)
        self.temp=self.temp.replace('[BOUND]','0 0 10 20')
        (Path.cwd()/'pes.cdxml').write_text(self.temp)
        # print(self.temp)

    