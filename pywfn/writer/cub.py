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

from pywfn.base import Mol
from pywfn import maths,base,config
from pywfn.utils import printer
from pywfn.spaceprop import wfnfunc,density
from pywfn.data.elements import elements
from pywfn import utils

class CubWriter:
    def __init__(self) -> None:
        """cube文件导出器"""
        self.title0='genetrate by pywfn'
        self.title1=time.strftime('%Y-%m-%d %H:%M:%S')
        self.syms:list[str]=[]         # 原子符号
        self.xyzs:np.ndarray|None=None # 原子坐标

        self.obts:list[int]=[]         # 分子轨道
        self.pos0:np.ndarray|None=None # 起点坐标
        self.size:list[int]  =[0,0,0]  # 三个方向的格点数量
        self.step:list[float]=[0,0,0]  # 三个方向的格点步长
        self.vals:np.ndarray|None=None # 格点数值[轨道,格点]
    
    def from_mol(self,mol:Mol):
        """从分子对象中读取数据"""
        self.syms=mol.atoms.syms
        self.xyzs=mol.atoms.xyzs
        return self
        
    def write_grids(self):
        """生成格点信息"""
        npos=self.size[0]*self.size[1]*self.size[2]
        self.file.write(f'{self.title0}\n{self.title1} {npos*len(self.obts)}\n')
        nx,ny,nz=self.size
        assert self.pos0 is not None,"没有设置起点坐标"
        x0,y0,z0=self.pos0
        sx,sy,sz=self.step
        natm=len(self.syms)
        self.file.write(f'{-natm:>5}{x0:>12.6f}{y0:>12.6f}{z0:>12.6f}\n')
        self.file.write(f'{nx:>5}{sx:>12.6f}{0:>12.6f}{0:>12.6f}\n')
        self.file.write(f'{ny:>5}{0:>12.6f}{sy:>12.6f}{0:>12.6f}\n')
        self.file.write(f'{nz:>5}{0:>12.6f}{0:>12.6f}{sz:>12.6f}\n')
    
    def write_coord(self):
        assert len(self.syms)!=0,"没有原子坐标"
        assert self.xyzs is not None,"没有原子坐标"
        assert len(self.syms)==len(self.xyzs),"原子坐标和符号数量不一致"
        natm=len(self.syms)
        for i in range(natm):
            x,y,z=self.xyzs[i]
            sym=self.syms[i]
            atomic=elements[sym].atomic
            self.file.write(f'{atomic:>5}{atomic:12.6f}{x:12.6f}{y:12.6f}{z:12.6f}\n')

    def write_value(self):
        """
        写入格点数值
        """
        nobt=len(self.obts)
        npos=self.size[0]*self.size[1]*self.size[2]
        assert self.vals is not None,"没有波函数值"
        assert utils.chkArray(self.vals,[nobt,npos]),f"数组形状不正确：{self.vals.shape}"
        index=0
        nx,ny,nz=self.size
        
        for i,info in enumerate([nobt]+self.obts): # 写入轨道信息
            self.file.write(f'{info:>5}')
            if(i+1)%10==0:self.file.write('\n')
        if(i+1)%10!=0:self.file.write('\n')

        for i in range(npos): # 对每一个点进行循环
            for j in range(nobt): # 对每一个轨道进行循环
                v=self.vals[j,i]
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