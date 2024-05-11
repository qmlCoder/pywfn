"""
将分子对象保存为gif文件
"""
from pathlib import Path
import json

from pywfn.base import Mol
from pywfn import config
from pywfn.data import temps

class gjfWriter:
    def __init__(self,mol:Mol):
        self.mol:Mol=mol
        self.temp=temps.gjf
        self.chkPath='chk'
        self.charge=0
        self.multi=1
        
    def get_coordStr(self):
        """生成坐标的字符串形式"""
        coordStrs=[]
        for atom in self.mol.atoms:
            x,y,z=atom.coord
            s=atom.symbol
            coordStr=f'{s:>2}{x:>14.8f}{y:>14.8f}{z:>14.8f}'
            coordStrs.append(coordStr)
        return '\n'.join(coordStrs)
    
    def search(self):
        """搜索需要的各种属性"""
        if self.mol.charge is not None:
            self.charge=self.mol.charge
        if self.mol.spin is not None:
            self.multi=self.mol.spin

    def build(self)->str:
        self.search()
        self.temp=self.temp.replace('<COORD>',self.get_coordStr())
        self.temp=self.temp.replace('<CHARGE>',f'{self.charge}')
        self.temp=self.temp.replace('<MULTI>',f'{self.multi}')
        self.temp=self.temp.replace('<CHK>',self.chkPath)
        return self.temp

    def save(self,path:str):
        self.build()
        Path(path).write_text(self.temp)