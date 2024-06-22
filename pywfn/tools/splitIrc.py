"""
分割IRC文件，IRC文件没有每一步的结构优化，所以只需要提取坐标即可
"""
import re,os
from pathlib import Path

from pywfn.writer import GjfWriter
from pywfn.data.elements import elements
from pywfn.utils import printer
from pywfn.base import Mol
from pywfn.reader import LogReader,AnyReader
from dataclasses import dataclass

import numpy as np

@dataclass
class Block:
    idx:int
    num:int

@dataclass
class AtmCord:
    sym:str
    x:float
    y:float
    z:float

@dataclass
class MolCord:
    atms:list[AtmCord]

class Tool:
    def __init__(self,path:str) -> None:
        self.path=Path(path)
        self.mol=Mol(LogReader(path))
        self.dirName=self.path.parent
        name=self.path.name # 包含后缀的文件名
        stem=self.path.stem # 不包含后缀的文件名
        self.fold=self.dirName / f'IRCS_{stem}'
        if not Path(self.fold).exists(): #判断文件夹是否存在
            os.mkdir(self.fold)
        self.content=self.path.read_text(encoding='utf-8')
        self.lines=self.content.splitlines(keepends=False)
        if 'Standard orientation' in self.content:
            self.coordType='Standard orientation'
        else:
            self.coordType='Input orientation'
        self.molCords:list[MolCord]=[]
        self.xyzs=[]
        self.syms=[]
        self.tits=[]
        self.mathc_titles()
        self.match_coords()
    
    def mathc_titles(self):
        for i,line in enumerate(self.lines): # 循环每一行
            if self.coordType in line: # 找到行坐标
                self.tits.append(i)
    
    def match_coords(self):
        s=r" +\d+ +(\d+) +\d+ +(-?\d+.\d+) +(-?\d+.\d+) +(-?\d+.\d+)"
        for tit in self.tits:
            xyzs=[]
            syms=[]
            for i in range(tit+5,len(self.lines)):
                line=self.lines[i]
                find=re.search(s,line)
                if find is None:break
                idx,x,y,z=find.groups()
                x,y,z=float(x),float(y),float(z)
                sym=elements[int(idx)].symbol
                xyzs.append([x,y,z])
                syms.append(sym)
            xyzs=np.array(xyzs)
            self.xyzs.append(xyzs)
            self.syms.append(syms)
        
    def build_mol(self,i:int):
        reader=AnyReader()
        reader.coords=self.xyzs[i]
        reader.symbols=self.syms[i]
        mol=Mol(reader)
        return mol
            


    def split(self):
        """分割文件"""
        # 首先获取所有构象的坐标
        for i in range(len(self.syms)):
            path=os.path.join(self.fold/f'f{i+1:0>2}.gjf')
            mol=self.build_mol(i)
            writer=GjfWriter(mol)
            writer.chkPath=f'f{i+1:0>2}'
            writer.save(path)