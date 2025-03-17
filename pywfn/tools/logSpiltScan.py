"""
输入文件有很多结构，将每个结构分割来来
重点在于分割，首先根据指定的分割符将文本分割为很多块，然后在每一块提取所需内容
重点在于提取结构信息，其他信息是不需要的，因此扫描的时候不要加入太多信息
"""
import re,os
from typing import Any
from pathlib import Path
from pywfn.data import Elements
from pywfn.utils import printer
from pywfn.base import Mole
from pywfn.writer import GjfWriter
from pywfn.reader import AnyReader
import numpy as np

class Tool:
    def __init__(self,path:str) -> None:
        """输入log文件的路径"""
        self.path=Path(path)
        self.fold=self.path.parent / self.path.stem
        
        self.content=self.path.read_text(encoding='utf-8')
        self.lines=self.content.splitlines(keepends=False)

        if 'Standard orientation' in self.content: # 坐标类型，是否为标准朝向
            self.coordType='Standard orientation'
        else:
            self.coordType='Input orientation'
        self.elements=Elements()
        self.title=''

    def split_raw(self)->list[str]:
        """将原始的一大段文本分割"""
        mark=' Optimization completed.' # 以此为每个结构的分隔，因为是柔性扫描，因此每一步结束都有这么一个标志
        # blocks=self.content.split(mark)
        blocks=re.split(mark,self.content)
        return blocks
    
    def find_titles(self,blocks:list[str]):
        """找到坐标对应的标题行"""
        titleNums=[]
        for idx,block in enumerate(blocks): # 循环每一段
            titleNum=[idx,0]
            lines=block.splitlines(keepends=False)
            for i,line in enumerate(lines): # 循环每一行
                
                if self.coordType in line: # 找到行坐标
                    titleNum=[idx,i]
                    # titleNums.append((idx,i)) # 只要第一个,记录段落数和行数，为啥只要第一个？不应该要最后一个？
                    # break
            titleNums.append(titleNum)
        return titleNums
    
    def match_coords(self,blocks:list[str],titleNums:list[tuple[int,int]])->tuple[list[str],np.ndarray]:
        """返回包含全部结构坐标的三维数组"""
        s=r" +\d+ +(\d+) +\d+ +(-?\d+.\d+) +(-?\d+.\d+) +(-?\d+.\d+)"
        xyzl=[]
        syml=[]
        for idx,titleNum in titleNums:
            xyzs=[]
            syms=[]
            lines=blocks[idx].splitlines(keepends=False)
            for i in range(titleNum+5,len(lines)): # 从标题后的第五行开始
                line=lines[i]
                search=re.search(s,line)
                if search is not None:
                    idx,x,y,z=search.groups()
                    symbol=self.elements[int(idx)].symbol
                    syms.append(symbol)
                    xyzs.append([float(x),float(y),float(z)])
                else:
                    break
            xyzl.append(xyzs)
            syml.append(syms)
        return syml,np.array(xyzl)*1.889
            
    def split(self):
        """分割文件"""
        blocks=self.split_raw()
        # 首先获取所有构象的坐标
        
        titleNums=self.find_titles(blocks)
        syml,xyzl=self.match_coords(blocks,titleNums)
        return syml,xyzl

    def save(self):
        """
        保存文件，可以保存到单个文件或多个文件
        """
        syml,xyzl=self.split()
        gjfStrs=[]
        for i in range(len(syml)-1):
            # 创建文件夹
            mol=Mole(AnyReader('',{
                'symbols':syml[i],
                'coords':xyzl[i],
                'charge':0,
                'spin':1,
            }))
            writer=GjfWriter().fromMol(mol)
            writer.title=self.title
            gjfStr=writer.build()
            gjfStrs.append(gjfStr)

        if not Path(self.fold).exists():os.mkdir(self.fold) # 创建文件夹
        for i,gjfStr in enumerate(gjfStrs):
            (self.fold/f's{i+1:0>2}.gjf').write_text(gjfStr,encoding='utf-8')
        
        printer.print("文件导出成功!!")

    def __call__(self) -> Any:
        self.blocks=self.split_raw()
