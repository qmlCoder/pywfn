"""
需要多个分子才能计算的属性
"""

from functools import cached_property,lru_cache
import numpy as np

from pywfn import base
from pywfn.atomprop import atomCharge
from pywfn.molprop import MolProp
from pywfn.utils import printer

class Fukui:
    """
    福井函数，以及依赖于福井函数计算的一些性质
    福井函数值针对原子的函数
    """
    def __init__(self,molN:"base.Mol",molM:"base.Mol",molP:"base.Mol") -> None:
        """
        需要三个分子对象
        molN:负电性分子
        molM:中性分子
        molP:正电性分子
        """
        
        if len(set([len(molN.atoms),len(molM.atoms),len(molP.atoms)]))!=1:
            printer.warn('三个分子的原子数量不相等')
        calerN=atomCharge.Calculator(molN)
        calerM=atomCharge.Calculator(molM)
        calerP=atomCharge.Calculator(molP)
        chargesN=calerN.calculate(atoms=None)
        chargesM=calerM.calculate(atoms=None)
        chargesP=calerP.calculate(atoms=None)
        self.fns=np.array(chargesM)-np.array(chargesN)
        self.fps=np.array(chargesP)-np.array(chargesM)
        self.mol=molM
        self.molProp=MolProp(molM)
    
    @lru_cache
    def fn(self,idx:int):
        """获取某个原子的福井函数+"""
        return self.fns[idx-1]
    
    @lru_cache
    def fp(self,idx:int):
        """获取某个原子的福井函数-"""
        return self.fps[idx-1]
    
    @lru_cache
    def cs(self,idx:int):
        """局部亲核指标"""
        cs=self.molProp.cs
        return cs*self.fn(idx),cs*self.fp(idx)
    
    @lru_cache
    def ei(self,idx:int):
        """局部亲电指标"""
        ei=self.molProp.ei
        return ei*self.fn(idx),ei*self.fp(idx)