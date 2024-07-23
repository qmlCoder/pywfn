"""
将分子对象保存为gif文件
"""
from pathlib import Path
import json

from pywfn.base import Mol
from pywfn import config
from pywfn.data import temps

class GjfWriter:
    def __init__(self,mol:Mol):
        self.mol:Mol=mol
        self.temp=temps.get('gjf')
        self.CHK='chk'
        self.CHARGE=0
        self.SPIN=1
        self.TITLE=''
        self.COORD=''
        
    def get_coordStr(self):
        """生成坐标的字符串形式"""
        coordStrs=[]
        for atom in self.mol.atoms:
            x,y,z=atom.coord/config.BOHR_RADIUS
            s=atom.symbol
            coordStr=f'{s:>2}{x:>14.8f}{y:>14.8f}{z:>14.8f}'
            coordStrs.append(coordStr)
        return '\n'.join(coordStrs)
    
    def search(self):
        """搜索需要的各种属性"""
        self.COORD=self.get_coordStr()
        try:
            self.CHARGE=self.mol.reader.get_charge()
        except:
            pass

        try:
            self.SPIN=self.mol.reader.get_spin()
        except:
            pass

    def build(self)->str:
        self.search()
        self.temp=self.temp.replace('<COORD>',self.COORD)
        self.temp=self.temp.replace('<CHARGE>',f'{self.CHARGE}')
        self.temp=self.temp.replace('<SPIN>',f'{self.SPIN}')
        self.temp=self.temp.replace('<CHK>',self.CHK)
        self.temp=self.temp.replace('<TITLE>',self.TITLE)
        return self.temp

    def save(self,path:str):
        self.build()
        Path(path).write_text(self.temp)
        print(f'文件导出至{path}')
    
    def onShell(self):
        path=Path(self.mol.reader.path)
        path=(path.parent/f'{path.stem}.gjf')
        self.save(f'{path}')
