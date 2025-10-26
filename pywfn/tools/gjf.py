"""
每个函数都是对gjf文件的修改或生成新的gjf，并返回新的gjf的内容
"""


from pathlib import Path
import numpy as np

from pywfn.base.mole import Mole
from pywfn.reader.gjf import GjfReader
from pywfn.writer.gjf import GjfWriter
from pywfn.editor import Editor
from pywfn.data import temps
from pywfn.data.elements import elements

class Tool:
    def __init__(self) -> None:
        pass

    def join(self,paths:list[str])->str:
        """将多个gjf文件拼接到一起"""
        texts=[]
        for path in paths:
            if Path(path).suffix!='.gjf':continue
            text=Path(path).read_text()
            texts.append(text)
        text:str='--Link1--\n\n'.join(texts)
        return text
    
    def ringBq(self,path:str,rings:list[list[int]],title:str='b3lyp/6-31g(d) NMR'):
        """在gjf的环中心添加Bq原子，方便计算NICS"""
        mol=Mole(GjfReader(path))
        for ring in rings:
            idxs=[atm-1 for atm in ring]
            x,y,z=np.mean(mol.xyzs[idxs,:],axis=0)
            mol.geome.addAtom('Bq',np.array([x,y,z]))
        writer=GjfWriter().fromMol(mol)
        writer.title=title
        writer.chk=f'{Path(mol.reader.path).stem}_Bq.chk'
        text=writer.build()
        return text
    
    # def scan_bond(self,path:str,atm1:int,atm2:int,step:int,size:float):
    #     reader=GjfReader(path)
    #     mole=Mole(reader)
    #     editor=Editor(mole)
    #     writer=GjfWriter().fromMol(mole)
    #     gjfs=[]
    #     for i in range(step):
    #         rmol=editor.rotate_bond(atm1,atm2,i*size)
    #         writer=GjfWriter().fromMol(rmol)
    #         gjfs.append(writer.build())
    #     text='--Link1--\n\n'.join(gjfs)
    #     return text
    
    def addElec(self,path:str,delta:int):
        """分子加减电子，影响电荷，加一个电子电荷变为-1

        Args:
            delEle (int): 加减电子数
        """
        reader=GjfReader(path)
        mol=Mole(reader)
        sat1=delta>0 # 增加还是减少
        nela,nelb=mol.nele
        for i in range(abs(delta)):
            sat2=nela==nelb
            match (sat1,sat2):
                case (True,True):
                    nela+=1
                case (True,False):
                    nelb+=1
                case (False,True):
                    nelb-=1
                case (False,False):
                    nela-=1
        writer=GjfWriter().fromMol(mol)
        atomic=sum(mol.atoms.atomics)
        writer.charge=atomic-nela-nelb
        writer.spin=nela-nelb+1
        return writer.build()
    