"""
将分子导出为xyz文件
"""
from pywfn import base
from pathlib import Path
from pywfn.data import elements
from pywfn.utils import printer

class xyzWriter:
    def __init__(self,mol:"base.Mol") -> None:
        self.mol=mol
        self.resStr=''

    def write(self):
        self.resStr+=f'{len(self.mol.atoms)}\n'
        self.resStr+=f'generate by pywfn\n'
        for atom in self.mol.atoms:
            x,y,z=atom.coord
            sym=atom.symbol
            idx=elements[sym].charge
            self.resStr+=f' {idx:<14}{x:>14.8f}{y:>14.8f}{z:>14.8f}\n'

    def save(self):
        self.write()
        path=self.mol.reader.path
        path=Path(path)
        (path.parent/f'{path.stem}.xyz').write_text(self.resStr)
        printer.res('导出成功!😄')