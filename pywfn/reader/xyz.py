"""
本脚本用来读取.xyz文件
xyz文件只包含原子坐标,所以很好读
一个xyz文件可以包含多个分子的坐标，但只生成最后一个坐标的分子对象
"""
from typing import List
import re
import numpy as np
from pywfn import base
from pywfn import reader

class XyzReader(reader.Reader):
    def __init__(self,path) -> None:
        self.path=path
        with open(path,'r',encoding='utf-8') as f:
            self.content=f.read()
        self.logLines=self.content.splitlines(keepends=False)
        self.geoms=self.read_geoms()

    def get_atmSyms(self) -> List[str]:
        return [sym for (sym,x,y,z) in self.geoms[-1]]
    
    def get_atmXyzs(self) -> np.ndarray:
        xyzs=[(x,y,z) for (sym,x,y,z) in self.geoms[-1]]
        return np.array(xyzs)
    
    def read_geoms(self):
        """
        有三种情况,原子数量,标题,原子坐标
        """
        
        p1=r'^ +\d+$'
        # 不满足第一种和第二种情况就是第三种情况
        p3=r'^ +([A-Za-z])+ +(-?\d+.\d+) +(-?\d+.\d+) +(-?\d+.\d+)$'
        geoms:list[list[tuple[str,float,float,float]]]=[]
        for line in self.logLines:
            if s1:=re.search(p1,line):
                geoms.append([])
            elif s3:=re.search(p3,line):
                sym,x,y,z=s3.groups()
                x,y,z=[float(x),float(y),float(z)]
                geoms[-1].append((sym,x,y,z))
            else:
                continue
        return geoms
