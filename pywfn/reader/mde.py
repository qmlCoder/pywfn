"""
读取.molden文件的脚本
"""
from pywfn import reader
from pywfn.base.basis import BasisData,Basis
from pywfn.base.coefs import Coefs
from pywfn.base.geome import Geome
from pywfn.reader.utils import toCart
from pywfn.data import bastrans
import numpy as np
import re
from functools import lru_cache
from pywfn import core


sym2ang={'s':0,'p':1,'d':2,'f':3,'g':4,'h':5,'i':6}

class MdeReader(reader.Reader):
    def __init__(self, path:str):
        super().__init__(path)
        self.type='mde'
        self.reader=core.reader.MdeReader(path) # type: ignore
    
    # def get_nele(self):
    #     engs,occs,CM=self.read_coefs()
    #     row,col=CM.shape
    #     if row>=col: # 闭壳层
    #         nela=nelb=sum(occs)
    #     else:
    #         nela=sum(occs[:col//2])
    #         nleb=sum(occs[col//2:])
    #     return nela,nelb
    
    def get_geome(self)->"Geome":
        """获取分子几何信息"""
        
        geome_core=self.reader.get_geome()
        geome=Geome()
        geome.core=geome_core
        return geome
    
    def get_basis(self)->"Basis":
        basis_core=self.reader.get_basis()
        basis=Basis()
        basis.core=basis_core
        return basis
    
    
    def get_coefs(self)->"Coefs":
        coefs_core=self.reader.get_coefs()
        coefs=Coefs()
        coefs.core=coefs_core
        return coefs