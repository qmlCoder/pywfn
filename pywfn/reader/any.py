"""
不从文件中读取信息的读取器
"""
import numpy as np

from pywfn import reader
from pywfn.data.basis import Basis
class AnyReader(reader.Reader):
    def __init__(self,path:str='',props:dict={}):
        super().__init__(path)
        self.charge:int=0
        self.spin:int=1
        self.symbols:list[str]=[]
        self.coords:np.ndarray

        for key,val in props.items(): # 根据传入的字典可以覆盖属性
            setattr(self,key,val)
    
    def get_charge(self) -> int:
        return self.charge
    
    def get_spin(self) -> int:
        return self.spin

    def get_symbols(self) -> list[str]:
        return self.symbols
    
    def get_coords(self) -> np.ndarray:
        return self.coords