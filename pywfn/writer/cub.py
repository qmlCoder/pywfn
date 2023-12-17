"""
用来生成格点数据并保存cube文件
cube文件使用的是玻尔半径B
分子的坐标使用的是埃米A
A=B/1,889
B=A*1.889
计算的时候使用A，写文件的时候将坐标和格点转为B
"""

import numpy as np
import time
from tqdm import tqdm
from pathlib import Path

from pywfn import maths,base,config
from pywfn.utils import printer
from pywfn.writer import lutils


class cubWriter:
    def __init__(self,mol:base.Mol) -> None:
        self.mol=mol
        self.step:float=config.RENDER_CLOUD_STEP
        self.border:float=config.RENDER_CLOUD_BORDER
        self.atoms:list[int]=[atom.idx for atom in self.mol.atoms]
        self.floatNum=12
        self.obts=[]
    
    def init_file(self):
        title0='genetrate by pywfn'
        title1=time.strftime('%Y-%m-%d %H:%M:%S')
        path=Path(self.mol.reader.path)
        obts=','.join((f'{e+1}' for e in self.obts))
        self.filePath=(path.parent/f'{path.stem}_{obts}.cub')
        self.file=open(self.filePath,mode='w')
        self.file.write(f'{title0}\n{title1}\n')

    def get_gridPos(self):
        """生成格点数据"""
        atoms=[atom-1 for atom in self.atoms]
        p0=self.mol.coords[atoms,:].min(axis=0)
        p1=self.mol.coords[atoms,:].max(axis=0)
        bord=self.border
        
        (Nx,Ny,Nz),gridPos=maths.gridPos(p0-bord,p1+bord,self.step,getN=True) #计算波函数时的坐标还是要使用原子单位的坐标
        assert Nx*Ny*Nz==len(gridPos),"网格点数量不匹配"
        x0,y0,z0=[e*config.BOHR_RADIUS for e in p0-bord]
        step=self.step*config.BOHR_RADIUS
        self.file.write(f'{-len(self.mol.atoms):>5}{x0:>12.6f}{y0:>12.6f}{z0:>12.6f}\n')
        self.file.write(f'{Nx:>5}{step:>12.6f}{0:>12.6f}{0:>12.6f}\n')
        self.file.write(f'{Ny:>5}{0:>12.6f}{step:>12.6f}{0:>12.6f}\n')
        self.file.write(f'{Nz:>5}{0:>12.6f}{0:>12.6f}{step:>12.6f}\n')
        self.write_coord()
        return gridPos
    

    
    def write_coord(self):
        for atom in self.mol.atoms:
            x,y,z=atom.coord*config.BOHR_RADIUS
            self.file.write(f'{atom.atomic:>5}{atom.atomic:12.6f}{x:12.6f}{y:12.6f}{z:12.6f}\n')
    
    def render(self):
        obts=self.obts
        gridPos=self.get_gridPos()
        self.file.write(f'{len(obts):>5}'+''.join(f'{obt+1:>5}' for obt in obts)+'\n')
        for obt in obts:
            values=self.mol.get_cloud(obt,gridPos,self.atoms)
            for i in printer.track(range(0,len(values),6),'正在写入'):
                self.file.write(''.join([f'{v:>13.5E}' for v in values[i:i+6]])+'\n')
    
    def save(self):
        self.init_file()
        self.render()
        printer.res(f'导出文件至{self.filePath}')
        self.file.close()