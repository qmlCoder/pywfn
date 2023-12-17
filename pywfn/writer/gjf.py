"""
将分子对象保存为gif文件
"""
from pathlib import Path
from pywfn.base import Mol
from pywfn import config
from pywfn.data import temps

class gjfWriter:
    def __init__(self,mol:Mol):
        self.mol:Mol=mol

        self.template=temps.gjf
        self.resStr=''
        
    def get_coordStr(self):
        """生成坐标的字符串形式"""
        coordStrs=[]
        for atom in self.mol.atoms:
            x,y,z=atom.coord
            s=atom.symbol
            coordStr=f'{s:>2}{x:>14.8f}{y:>14.8f}{z:>14.8f}'
            coordStrs.append(coordStr)
        return '\n'.join(coordStrs)
    
    def save(self):
        path=Path(self.mol.reader.path)
        replaces=[
            ['<COORD>', self.get_coordStr()],
            ['<CHK>', str(path.stem)],
            ['<CHARGE>',f'{self.mol.charge}'],
            ['<MULTI>',f'{self.mol.spin}']
        ]
        content=self.template
        for k,v in replaces:
            content=content.replace(k,v)
        print(content)
        resFold=path.parent /f'gjfs'
        if not resFold.exists():resFold.mkdir()
        resPath=resFold/f'{path.stem}.gjf'
        resPath.write_text(content)