"""
将分子对象保存为gif文件
"""
from pathlib import Path
import numpy as np
from pywfn.base.mole import Mole
from pywfn import config
from pywfn.datas import temps
from pywfn import core

class GjfWriter:
    def __init__(self):
        self.core=core.writer.GjfWriter() # type: ignore

    def set_syms(self,syms:list[str]):
        self.core.set_syms(syms)
    
    def set_xyzs(self,xyzs:np.ndarray):
        self.core.set_xyzs(xyzs)

    def set_chk(self,chk:str):
        self.core.set_chk(chk)

    def set_job(self,job:str):
        self.core.set_job(job)
    
    def set_charge(self,charge:int):
        self.core.set_charge(charge)
    
    def set_multip(self,multip:int):
        self.core.set_multip(multip)
    
    def from_mole(self,mole:Mole):
        self.core.from_mole(mole.core)
    
    def build(self)->str:
        return self.core.build()

    def save(self,path:str):
        return self.core.save(path)