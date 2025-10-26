from pywfn.base import Mole
from typing import Literal
from pywfn import core
from pywfn.moleprop.orbital import Deco

Chrgs=Literal['mulliken','lowdin','space','hirshfeld']

class Calculator():
    def __init__(self,mole:"Mole"):
        self.logTip:str=''
        self.mole=mole
        self.form:str='charge' # 输出格式，电子数或电荷数 number|charge
        self.caler=core.atomprop.charge.Calculator(mole.mole) # type: ignore # 核心计算器

    def spin(self,ctype:str='mulliken'):
        """计算所有原子的自旋"""
        return self.caler.spin()
    
    def mulliken(self,form:str='charge'):
        """
        计算mulliken电荷
        """
        return self.caler.mulliken()
    
    def lowdin(self,form:str='charge'):
        """
        计算每个原子的lowdin电荷
        """
        return self.caler.lowdin()
    
    def hirshfeld(self,form:str='charge'):
        """计算原子的Hirshfeld电荷"""
        return self.caler.hirshfeld()

    def pocv(self,dirs:dict[int,list[float]],keep_other_atm:bool,keep_other_sym:bool,ctype:str):
        """计算投影分子轨道的电子数"""
        return self.caler.pocv(dirs,keep_other_atm,keep_other_sym,ctype) # type: ignore
    
    def deco(self,decos:dict[int,Deco],ctype:str='mulliken'):
        return self.caler.deco(decos,ctype)
    
    def pi_pocv(self,ctype:str='mulliken',atms:list[int]|None=None):
        """
        计算π电子数
        每个原子的方向为`法向量`的`方向电子数`即为`π电子数`
        """
        return self.caler.pi_pocv(ctype,atms)
    
    def pi_deco(self,ctype='mulliken'): # 使用轨道分解方法计算pi电子分布，可以包含D轨道
        return self.caler.pi_deco(ctype)
