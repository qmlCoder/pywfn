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
from pathlib import Path

from pywfn import maths,base,config
from pywfn.utils import printer
from pywfn.spaceProp import wfnfunc,density

class CubWriter:
    def __init__(self,mol:base.Mol) -> None:
        self.mol=mol
        self.step:float=config.RENDER_CLOUD_STEP
        self.border:float=config.RENDER_CLOUD_BORDER
        self.atms:list[int]=self.mol.atoms.atms
        self.direct:np.ndarray=None
        self.floatNum=12
        na,nb=self.mol.eleNum
        self.obts=[na-1]
        if self.mol.open:
            nbas=self.mol.CM.shape[0]
            self.obts.append(nbas+nb-1)
        self.wfnCaler=wfnfunc.Calculator(mol)
        self.denCaler=density.Calculator(mol)
        self.CM=self.mol.CM
        self.ctype='wfn' # 计算的类型，波函数、电子密度
    
    def init_file(self,path:str):
        self.title0='genetrate by pywfn'
        self.title1=time.strftime('%Y-%m-%d %H:%M:%S')
        # obts=','.join((f'{e+1}' for e in self.obts))
        self.filePath=path
        self.file=open(self.filePath,mode='w')
        
    def get_gridPos(self):
        """生成格点数据"""
        self.mol.bohr=False
        idxs=[atm-1 for atm in self.atms]
        p0=self.mol.coords[idxs,:].min(axis=0)
        p1=self.mol.coords[idxs,:].max(axis=0)
        bord=self.border
        
        (Nx,Ny,Nz),gridPos=maths.cubeGrid(p0,p1,self.step,bord=bord) #计算波函数时的坐标还是要使用原子单位的坐标
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

    def calc_values(self):
        """
        计算数值，计算和保存逻辑分开
        """
    
    def write_value(self):
        """
        计算波函数值并写入文件
        """
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
            self.wfnCaler.set_grid(gridPos)
            self.wfnCaler.atms=self.atms
            allVals=[self.wfnCaler.obtWfn(obt) for obt in obts]
        elif self.ctype=='den':
            self.denCaler.set_grid(gridPos)
            allVals=[self.denCaler.molDens_lib()]
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
    
    def save(self,path:str):
        self.init_file(path)
        self.write_value()
        printer.res(f'导出文件至{self.filePath}')
        self.file.close()

    def onShell(self):
        name='开窍层' if self.mol.open else '闭壳层'
        na,nb=self.mol.eleNum
        print(f'该分子为{name},alpha和beta电子数量分别为{na},{nb}')
        from pywfn.utils import parse_intList,parse_obtList
        atms=input('输入需渲染的原子(默认所有原子): ')
        if atms!='':self.atoms=parse_intList(atms,start=1)
        obts=input('输入需渲染的轨道(默认HOMO轨道),例如a10,b10: ')
        nbas=self.mol.CM.shape[0]
        if obts!='':self.obts=parse_obtList(obts,nbas)
        path=self.mol.reader.path
        self.save(f'{path}.cub')