"""
分割IRC文件，IRC文件没有每一步的结构优化，所以只需要提取坐标即可
"""
import re,os
from pathlib import Path

from pywfn.tools.fileCreater import Tool as  FileCreater
from pywfn.data import Elements
from pywfn.utils import printer

class Tool:
    def __init__(self,path:str) -> None:
        self.path=Path(path)
        self.dirName=self.path.parent
        self.baseName=self.path.name # 包含后缀的文件名
        self.fileName=self.path.stem # 不包含后缀的文件名
        folder=self.dirName / self.fileName
        if not Path(folder).exists(): #判断文件夹是否存在
            os.mkdir(folder)
        self.content=self.path.read_text(encoding='utf-8')
        self.lines=self.content.splitlines(keepends=False)
        if 'Standard orientation' in self.content:
            self.coordType='Standard orientation'
        else:
            self.coordType='Input orientation'
        self.elements=Elements()
    
    def find_titles(self):
        titleNums=[]
        idx=0
        for i,line in enumerate(self.lines): # 循环每一行
            if self.coordType in line: # 找到行坐标
                titleNums.append((idx,i)) # 只要第一个,记录段落数和行数
                idx+=1
        return titleNums
    
    def match_coords(self,titleNums:list[tuple[int]]):
        s=r" +\d+ +(\d+) +\d+ +(-?\d+.\d+) +(-?\d+.\d+) +(-?\d+.\d+)"
        coords=[]
        for idx,titleNum in titleNums:
            coord=[]
            for i in range(titleNum+5,len(self.lines)):
                line=self.lines[i]
                if re.search(s,line) is not None:
                    idx,x,y,z=re.search(s,line).groups()
                    symbol=self.elements.get_element_by_idx(int(idx)).symbol
                    coord.append([symbol,x,y,z])
                else:
                    break
            coords.append(coord)
        return coords

    def split(self):
        """分割文件"""
        # 首先获取所有构象的坐标
        titleNums=self.find_titles()
        coords=self.match_coords(titleNums)
        for i,coord in printer.track(enumerate(coords)):
            path=os.path.join(self.dirName,self.fileName,f'f{i+1:0>2}.gjf')
            fileCreater=FileCreater(path=path)
            fileCreater.set_coord(coord)
            fileCreater.set_chk(f'{i+1}')
            fileCreater.save()