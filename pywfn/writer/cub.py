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
from pywfn.spaceProp import wfnfunc,density

class CubWriter:
    def __init__(self,mol:base.Mol) -> None:
        self.mol=mol
        self.step:float=config.RENDER_CLOUD_STEP
        self.border:float=config.RENDER_CLOUD_BORDER
        self.atoms:list[int]=self.mol.atoms.indexs
        self.direct:np.ndarray=None
        self.floatNum=12
        self.obts=self.mol.O_obts
        self.wfnCaler=wfnfunc.Calculator(mol)
        self.denCaler=density.Calculator(mol)
        self.CM=self.mol.CM
        self.ctype='wfn' # 计算的类型，波函数、电子密度
    
    def init_file(self,name:str):
        self.title0='genetrate by pywfn'
        self.title1=time.strftime('%Y-%m-%d %H:%M:%S')
        path=Path(self.mol.reader.path)
        # obts=','.join((f'{e+1}' for e in self.obts))
        self.filePath=(path.parent/f'{path.stem}_{name}.cub')
        self.file=open(self.filePath,mode='w')
        
    def get_gridPos(self):
        """生成格点数据"""
        self.mol.bohr=False
        atoms=[atom-1 for atom in self.atoms]
        p0=self.mol.coords[atoms,:].min(axis=0)
        p1=self.mol.coords[atoms,:].max(axis=0)
        bord=self.border
        
        (Nx,Ny,Nz),gridPos=maths.gridPos(p0-bord,p1+bord,self.step) #计算波函数时的坐标还是要使用原子单位的坐标
        self.gridSize=(Nx,Ny,Nz)
        assert Nx*Ny*Nz==len(gridPos),"网格点数量不匹配"
        self.file.write(f'{self.title0}\n{self.title1} {len(gridPos)*len(self.obts)}\n')
        x0,y0,z0=[e for e in p0-bord]
        step=self.step
        natm=len(self.mol.atoms)
        self.file.write(f'{-natm:>5}{x0:>12.6f}{y0:>12.6f}{z0:>12.6f}\n')
        self.file.write(f'{Nx:>5}{step:>12.6f}{0:>12.6f}{0:>12.6f}\n')
        self.file.write(f'{Ny:>5}{0:>12.6f}{step:>12.6f}{0:>12.6f}\n')
        self.file.write(f'{Nz:>5}{0:>12.6f}{0:>12.6f}{step:>12.6f}\n')
        self.write_coord()
        return gridPos
    
    def write_coord(self):
        for atom in self.mol.atoms:
            x,y,z=atom.coord
            self.file.write(f'{atom.atomic:>5}{atom.atomic:12.6f}{x:12.6f}{y:12.6f}{z:12.6f}\n')
    
    def write_value(self):
        obts=self.obts
        gridPos=self.get_gridPos()
        Nx,Ny,Nz=self.gridSize
        if self.ctype=='wfn':
            obtInfo=[len(obts)]+[o+1 for o in obts]
        elif self.ctype=='den':
            obtInfo=[len(obts)]+[0]
        
        for i,info in enumerate(obtInfo): # 写入轨道信息
            self.file.write(f'{info:>5}')
            if(i+1)%10==0:self.file.write('\n')
        if(i+1)%10!=0:self.file.write('\n')
        if self.ctype=='wfn':
            allVals=[self.wfnCaler.obtWfn(obt,gridPos,atms=self.atoms,coefs=self.CM[:,obt]) for obt in obts]
        elif self.ctype=='den':
            allVals=[self.denCaler.molDens_obt(gridPos,atms=self.atoms,CM=self.CM,obts=self.obts)]
        lenVals=len(gridPos)
        index=0
        for i in printer.track(range(lenVals)): # 对每一个点进行循环
            for j in range(len(allVals)): # 对每一个轨道进行循环
                v=allVals[j][i]
                # if v==0:v=1e-8
                self.file.write(f'{v:13.5E}')
                index+=1
                if index==Nz*len(obts):
                    self.file.write('\n')
                    index=0
                    continue
                if index%6==0:self.file.write('\n')
    
    def save(self,name):
        self.init_file(name)
        self.write_value()
        printer.res(f'导出文件至{self.filePath}')
        self.file.close()

    def onShell(self):
        from pywfn.utils import parse_intList
        atms=input('输入要导出的原子[*]')
        if atms!='':self.atoms=parse_intList(atms)
        step=input(f'输入格点步长[{self.step:>.2f}]')
        if step!='':step=float(step)
        self.save(self.mol.reader.fname)