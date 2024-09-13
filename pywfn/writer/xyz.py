"""
将分子导出为xyz文件
"""
from pywfn import base
from pathlib import Path
from pywfn.data.elements import elements
from pywfn.utils import printer
from pywfn.data import temps
from pywfn import config

class XyzWriter:
    def __init__(self,mol:"base.Mol") -> None:
        self.mol=mol
        self.title='generate by pywfn'
        self.temp=temps.xyz
        self.atomForm='sym' # 打印原子的类型，元素符号[sym]或核电荷数[idx]

    def build(self):
        natm=len(self.mol.atoms)
        self.temp=self.temp.replace('<NATM>',f'{natm}')
        self.temp=self.temp.replace('<TITLE>',self.title)
        
        coordStrs=[]
        for atom in self.mol.atoms:
            x,y,z=atom.coord/config.BOHR_RADIUS
            sym=atom.symbol
            
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

    def onShell(self):
        path=Path(self.mol.reader.path)
        path=(path.parent/f'{path.stem}.xyz')
        self.save(path)