"""
本脚本用来读取.xyz文件
xyz文件只包含原子坐标,所以很好读
一个xyz文件可以包含多个分子的坐标,所以要生成多个分子对象
"""
from typing import List
import re

from pywfn import base
from pywfn import reader

class XyzReader(reader.Reader):
    def __init__(self,path) -> None:
        self.path=path
        with open(path,'r',encoding='utf-8') as f:
            self.content=f.read()
        self.logLines=self.content.splitlines(keepends=False)
        self.mols:List["base.Mole"]=[]
        self.read_coord()
    
    def read_coord(self):
        """
        有三种情况,原子数量,标题,原子坐标
        """
        
        s1='^ +\d+$'
        # 不满足第一种和第二种情况就是第三种情况
        s3='^ +([A-Za-z])+ +(-?\d+.\d+) +(-?\d+.\d+) +(-?\d+.\d+)$'

        for line in self.logLines:
            if re.search(s1,line) is not None:
                self.mols.append(base.Mole())
            elif re.search(s3,line) is not None:
                symbol,x,y,z=re.search(s3,line).groups()
                self.mols[-1].add_atom(symbol,[float(x),float(y),float(z)])
            else:
                continue
