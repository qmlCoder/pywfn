"""
在gjf的环中心添加Bq原子，方便计算NICS
"""
from pywfn.base import Mole
from pywfn.cli import Shell
from pywfn.writer import GjfWriter

import numpy as np
from pathlib import Path

class Tool:
    def __init__(self,mol:Mole) -> None:
        self.mol=mol
        self.rings=[]
        self.writer=GjfWriter().fromMol(mol)
        self.gjftxt=''
        self.title='b3lyp/6-31g(d) NMR' # 计算NICS只要个单点就行了

    def build(self):
        """构建"""
        molCoords=self.mol.coords.copy()
        for ring in self.rings:
            idxs=[atm-1 for atm in ring]
            cords=molCoords[idxs,:]
            coord=np.mean(cords,axis=0)
            self.mol._atoms.add('Bq',coord)
        self.writer.title=self.title
        self.writer.chk=f'{Path(self.mol.reader.path).stem}_Bq.chk'
        return self.writer.build()

    def save(self,path:str=''):
        if not path:
            fold=Path(self.mol.reader.path).parent
            stem=Path(self.mol.reader.path).stem
            path=f'{fold}/{stem}_Bq.gjf'
        self.gjftxt=self.build()
        Path(path).write_text(self.gjftxt)
        print(f'文件导出至{path}')