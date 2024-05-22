"""
定义mol文件的写出器
"""
from pywfn.base import Mol
from pywfn.data import temps
from pathlib import Path

class MolWriter:
    def __init__(self,mol:Mol):
        self.mol = mol
        self.temp=temps.mol
        self.title:str=None

    def atomBlock(self)->str:
        atomLines=[]
        for atom in self.mol.atoms:
            x,y,z=atom.coord/1.889
            s=atom.symbol
            line=f'{x:>10.4f}{y:>10.4f}{z:>10.4f}{s:>2}   0  0  0  0  0  0  0  0  0  0  0  0'
            atomLines.append(line)
        return '\n'.join(atomLines)
    
    def bondBlock(self)->str:
        bondLines=[]
        for bond in self.mol.bonds:
            a1,a2=bond.ats
            line=f'{a1:>3d}{a2:>3d}{bond.btype:>3d}  0  0  0  0'
            bondLines.append(line)
        return '\n'.join(bondLines)
    
    def build(self)->str:
        natom=self.mol.atoms.num
        nbond=self.mol.bonds.num
        if self.title is None:self.title=self.mol.formula
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
    
    def onShell(self):
        path=Path(self.mol.reader.path)
        path=(path.parent/f'{path.stem}.mol')
        self.save(path)
