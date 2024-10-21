"""
给gjf文件添加或减少一个电子
"""

from pywfn.base import Mol
from pywfn.writer import GjfWriter
from pathlib import Path
class Tool():
    def __init__(self,mol:Mol) -> None:
        self.mol=mol
    
    def addEle(self,delEle:int):
        """分子加减电子，影响电荷，加一个电子电荷变为-1

        Args:
            delEle (int): 加减电子数
        """
        na,nb=self.mol.eleNum
        # print('alpha和beta电子数分别为',na,nb)
        if na==nb:
            na+=1
        elif na>nb:
            nb+=1
        # print('alpha和beta电子数分别为',na,nb)
        charge=sum(self.mol.atoms.atomics)-na-nb
        spin=na-nb+1
        self.mol.props.set('charge',charge)
        self.mol.props.set('spin',spin)
    
    def save(self,path:str):
        writer=GjfWriter(self.mol)
        writer.CHK=Path(path).stem
        writer.save(path)