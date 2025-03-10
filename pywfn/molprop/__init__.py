"""
计算分子属性
"""
from functools import cached_property
from pywfn import base
from pywfn.utils import printer

# 记录TCE的HOMO能量
TCEs={
    '6-31G*':-0.33520
}

class MolProp:
    def __init__(self,mol:"base.Mol") -> None:
        self.mol=mol
    
    @cached_property
    def hlIdx(self):
        """homo和lomo的索引"""
        oets=self.mol.obtOccs
        for i,e in enumerate(oets):
            if e:continue
            return i-1,i
        raise Exception("没有找到homo和lomo")

    @cached_property
    def homo(self):
        """HOMO轨道的能量"""
        
        h,l=self.hlIdx
        return self.mol.obtEngs[h]
    
    @cached_property
    def lomo(self):
        """HOMO轨道的能量"""
        h,l=self.hlIdx
        return self.mol.obtEngs[l]
    
    @cached_property
    def ecp(self):
        """化学势"""
        return (self.homo+self.lomo)/2
    
    @cached_property
    def en(self):
        """电负性"""
        return -self.ecp

    @cached_property
    def ch(self):
        """化学硬度"""
        return self.lomo-self.homo
    
    @cached_property
    def cs(self):
        """化学软度"""
        return 1/self.ch
    
    @cached_property
    def ei(self):
        """亲电指数"""
        return self.ecp**2/self.ch
    
    @cached_property
    def N(self):
        """亲核反应指标"""
        basi=self.mol.basis.name

        return self.homo-TCEs[basi]

    def props(self):
        """以字典的形式返回分子属性"""
        props = [
            ['HOMO',self.homo],
            ['LOMO',self.lomo],
            ['化学势',self.ecp],
            ['电负性',self.en],
            ['化学硬度',self.ch],
            ['化学软度',self.cs],
            ['亲点指数',self.ei],
            ['亲和反应指标',self.N]
        ]
        for k,v in props:
            printer.info(f'{k:{chr(12288)}<6}\t{v:>10.6f}')