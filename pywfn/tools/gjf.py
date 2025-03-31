from pathlib import Path
import numpy as np

from pywfn.base import Mole
from pywfn.writer import GjfWriter

class Tool:
    def __init__(self) -> None:
        pass

    def join(self,paths:list[str],path:str='')->str:
        """将多个gjf文件拼接到一起"""
        texts=[]
        for path in paths:
            if Path(path).suffix!='.gjf':continue
            text=Path(path).read_text()
            texts.append(text)
        text:str='--Link1--\n\n'.join(texts)
        if path=='':
            return text
        else:
            Path(path).write_text('--Link1--\n\n'.join(texts),encoding='utf-8')
        return text
    
    def ringBq(self,mol:Mole,rings:list[list[int]],title:str='b3lyp/6-31g(d) NMR',path:str=''):
        """在gjf的环中心添加Bq原子，方便计算NICS"""
        writer=GjfWriter().fromMol(mol)
        molCoords=mol.coords.copy()
        for ring in rings:
            idxs=[atm-1 for atm in ring]
            x,y,z=np.mean(molCoords[idxs,:],axis=0)
            mol.geome.addAtom('Bq',np.array([x,y,z]))
        writer.title=title
        writer.chk=f'{Path(mol.reader.path).stem}_Bq.chk'
        text=writer.build()
        if path=='':
            return text
        else:
            Path(path).write_text(text,encoding='utf-8')
        return text
    