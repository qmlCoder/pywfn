"""
将分子导出为xyz文件
"""
from pywfn import base
from pywfn.base import Mol
from pathlib import Path
from pywfn.data.elements import elements
from pywfn.utils import printer
from pywfn.data import temps
from pywfn import config
import numpy as np

class XyzWriter:
    def __init__(self) -> None:
        self.title='generate by pywfn'
        self.temp=temps.get('xyz')
        self.atomForm='sym' # 打印原子的类型，元素符号[sym]或核电荷数[idx]
        self.syms:list[str]=[]    # 元素符号
        self.xyzs:np.ndarray|None=None # 原子坐标

    def fromMol(self,mol:Mol):
        self.syms=mol.atoms.syms
        self.xyzs=mol.atoms.xyzs
        return self

    def build(self):
        assert self.xyzs is not None,"没有提供坐标"
        natm=len(self.syms)
        self.temp=self.temp.replace('<NATM>',f'{natm}')
        self.temp=self.temp.replace('<TITLE>',self.title)
        
        coordStrs=[]
        for sym,xyz in zip(self.syms,self.xyzs):
            x,y,z=xyz/config.BOHR_RADIUS
            
            if self.atomForm=='sym':
                coordStrs.append(f' {sym:<14}{x:>14.8f}{y:>14.8f}{z:>14.8f}')
            elif self.atomForm=='idx':
                idx=elements[sym].charge
                coordStrs.append(f' {idx:<14}{x:>14.8f}{y:>14.8f}{z:>14.8f}')
            else:
                raise ValueError('atomForm must be sym or idx')
        self.temp=self.temp.replace('<COORD>','\n'.join(coordStrs))

    def save(self,path:str):
        self.build()
        Path(path).write_text(self.temp)
        printer.res(f'{path} 导出成功!😄')