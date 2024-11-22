"""
定义mol文件的写出器
"""
from pywfn.base import Mol
from pywfn.data import temps
from pathlib import Path
import numpy as np

class MolWriter:
    def __init__(self):
        self.temp=temps.mol
        self.title:str='title'
        self.syms:list[str]=[] #原子符号列表
        self.xyzs:np.ndarray|None=None #坐标数组
        self.bonds:list[tuple[int,int,int]]=[] #键列表 (原子1,原子2,键类型) 1：单键，2：双键，3：方向键

    def fromMol(self,mol:Mol):
        self.title=mol.formula
        self.syms=mol.atoms.syms
        self.xyzs=mol.atoms.xyzs
        self.bonds=[(bond.a1.idx,bond.a2.idx,bond.btype) for bond in mol.bonds]
        # self.bonds=[(1,2,3)]
        return self

    def atomBlock(self)->str:
        assert self.xyzs is not None,"没有提供坐标"
        atomLines=[]
        for sym,xyz in zip(self.syms,self.xyzs):
            x,y,z=xyz/1.889
            line=f'{x:>10.4f}{y:>10.4f}{z:>10.4f}{sym:>2}   0  0  0  0  0  0  0  0  0  0  0  0'
            atomLines.append(line)
        return '\n'.join(atomLines)
    
    def bondBlock(self)->str:
        bondLines=[]
        for a1,a2,btype in self.bonds:
            line=f'{a1:>3d}{a2:>3d}{btype:>3d}  0  0  0  0'
            bondLines.append(line)
        return '\n'.join(bondLines)
    
    def build(self)->str:
        natom=len(self.syms)
        nbond=len(self.bonds)
        self.temp=self.temp.replace('[TITLE]',self.title)
        self.temp=self.temp.replace('[NATOM]',f'{natom:>3d}')
        self.temp=self.temp.replace('[NBOND]',f'{nbond:>3d}')
        self.temp=self.temp.replace('[ATOMBLOCK]',self.atomBlock())
        self.temp=self.temp.replace('[BONDBLOCK]',self.bondBlock())
        return self.temp
    
    def save(self,path:str):
        """
        保存mol文件到指定路径
        path:str 文件路径
        """
        self.build()
        Path(path).write_text(self.temp)