"""
用来生成格点数据并保存cube文件
cube文件使用的是玻尔半径B
分子的坐标使用的是埃米A
A=B/1.889
B=A*1.889
计算的时候使用A，写文件的时候将坐标和格点转为B
"""

import numpy as np
import time
from pathlib import Path

from pywfn import maths,base,config
from pywfn.utils import printer
from pywfn.spaceProp import wfnfunc,density
from pywfn.data.elements import elements

class CubWriter:
    def __init__(self,syms:list[str],xyzs:np.ndarray,obts:list[int],pos0:list[float],nums:list[int],step:list[float],vals:list[np.ndarray]) -> None:
        """cube文件导出器

        Args:
            syms (list[str]): 原子符号
            xyzs (np.ndarray): 原子坐标
            obts (list[int]): 原子轨道
            pos0 (list[float]): 起点坐标
            nums (list[int]): 格点数量
            step (list[float]): 格点步长
            vals (list[np.ndarray]): 格点数量
        """
        self.title0='genetrate by pywfn'
        self.title1=time.strftime('%Y-%m-%d %H:%M:%S')
        self.syms=syms
        self.xyzs=xyzs
        self.obts=obts
        self.pos0=pos0
        self.nums=nums
        self.step=step
        self.vals=vals
        self.npos=self.nums[0]*self.nums[1]*self.nums[2]
        
    def write_grids(self):
        """生成格点数据"""
        
        self.file.write(f'{self.title0}\n{self.title1} {self.npos*len(self.obts)}\n')
        nx,ny,nz=self.nums
        x0,y0,z0=self.pos0
        sx,sy,sz=self.step
        natm=len(self.syms)
        self.file.write(f'{-natm:>5}{x0:>12.6f}{y0:>12.6f}{z0:>12.6f}\n')
        self.file.write(f'{nx:>5}{sx:>12.6f}{0:>12.6f}{0:>12.6f}\n')
        self.file.write(f'{ny:>5}{0:>12.6f}{sy:>12.6f}{0:>12.6f}\n')
        self.file.write(f'{nz:>5}{0:>12.6f}{0:>12.6f}{sz:>12.6f}\n')
    
    def write_coord(self):
        natm=len(self.syms)
        for i in range(natm):
            x,y,z=self.xyzs[i]
            sym=self.syms[i]
            atomic=elements[sym].atomic
            self.file.write(f'{atomic:>5}{atomic:12.6f}{x:12.6f}{y:12.6f}{z:12.6f}\n')

    
    def write_value(self):
        """
        计算波函数值并写入文件
        """
        nobt=len(self.obts)
        npos=self.npos
        index=0
        nx,ny,nz=self.nums

        
        for i,info in enumerate([nobt]+self.obts): # 写入轨道信息
            self.file.write(f'{info:>5}')
            if(i+1)%10==0:self.file.write('\n')
        if(i+1)%10!=0:self.file.write('\n')

        for i in range(npos): # 对每一个点进行循环
            for j in range(nobt): # 对每一个轨道进行循环
                v=self.vals[j][i]
                # if v==0:v=1e-8
                self.file.write(f'{v:13.5E}')
                index+=1
                if index==nz*nobt:
                    self.file.write('\n')
                    index=0
                    continue
                if index%6==0:self.file.write('\n')
    
    def save(self,path:str):
        self.file=open(path,mode='w')
        self.write_grids()
        self.write_coord()
        self.write_value()
        printer.res(f'导出文件至{path}')
        self.file.close()