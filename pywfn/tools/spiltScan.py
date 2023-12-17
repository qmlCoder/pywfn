"""
输入文件有很多结构，将每个结构分割来来
重点在于分割，首先根据指定的分割符将文本分割为很多块，然后在每一块提取所需内容
"""
import re,os
from typing import Any
from pathlib import Path

from pywfn.tools import writer
from pywfn.data import Elements
from pywfn.utils import printer

class Tool:
    def __init__(self,path:str) -> None:
        """输入log文件的路径"""
        self.path=Path(path)
        self.dirName=self.path.parent
        self.baseName=self.path.name # 包含后缀的文件名
        self.fileName=self.path.stem # 不包含后缀的文件名
        folder=self.dirName / self.fileName
        if not Path(folder).exists(): #判断文件夹是否存在
            os.mkdir(folder) # 创建文件夹
        self.content=self.path.read_text(encoding='utf-8')
        self.lines=self.content.splitlines(keepends=False)

        if 'Standard orientation' in self.content: # 坐标类型，是否为标准朝向
            self.coordType='Standard orientation'
        else:
            self.coordType='Input orientation'
        self.elements=Elements()

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
    
    def match_coords(self,blocks:list[str],titleNums:list[tuple[int]])->list[list[float]]:
        """返回包含全部结构坐标的三维数组"""
        s=r" +\d+ +(\d+) +\d+ +(-?\d+.\d+) +(-?\d+.\d+) +(-?\d+.\d+)"
        coords=[]
        for idx,titleNum in titleNums:
            coord=[]
            lines=blocks[idx].splitlines(keepends=False)
            for i in range(titleNum+5,len(lines)): # 从标题后的第五行开始
                line=lines[i]
                if re.search(s,line) is not None:
                    idx,x,y,z=re.search(s,line).groups()
                    symbol=self.elements[int(idx)].symbol
                    coord.append([symbol,float(x),float(y),float(z)])
                else:
                    break
            coords.append(coord)
        return coords
            
    def split(self):
        """分割文件"""
        blocks=self.split_raw()
        # 首先获取所有构象的坐标

        titleNums=self.find_titles(blocks)
        coords=self.match_coords(blocks,titleNums)
        for i,coord in printer.track(enumerate(coords)):
            path=os.path.join(self.dirName,self.fileName,f'f{i+1:0>2}.gjf')
            writer.gjf(path,f'{i+1}',coord,0,1)

    def __call__(self) -> Any:
        self.blocks=self.split_raw()
