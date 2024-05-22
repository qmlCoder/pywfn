"""
将分子导出为xyz文件
"""
from pywfn import base
from pathlib import Path
from pywfn.data.elements import elements
from pywfn.utils import printer
from pywfn.data import temps

class XyzWriter:
    def __init__(self,mol:"base.Mol") -> None:
        self.mol=mol
        self.title='generate by pywfn'
        self.temp=temps.xyz

    def build(self):
        natm=len(self.mol.atoms)
        self.temp=self.temp.replace('<NATM>',f'{natm}')
        self.temp=self.temp.replace('<title>',self.title)
        
        coordStrs=[]
        for atom in self.mol.atoms:
            x,y,z=atom.coord
            sym=atom.symbol
            idx=elements[sym].charge
            coordStrs.append(f' {idx:<14}{x:>14.8f}{y:>14.8f}{z:>14.8f}')
        self.temp=self.temp.replace('<coord>','\n'.join(coordStrs))

    def save(self,path:str):
        self.build()
        Path(path).write_text(self.temp)
        printer.res(f'{path} 导出成功!😄')

    def onShell(self):
        path=Path(self.mol.reader.path)
        path=(path.parent/f'{path.stem}.xyz')
        self.save(path)