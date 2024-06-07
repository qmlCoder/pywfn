"""
实现单位转换
既然作为波函数分析的程序，针对的都是原子和分子，那么使用原子单位就比较合适。
使用原子单位 https://zhuanlan.zhihu.com/p/362046144
"""

# 定义支持的单位和类型，同一个类型的单位之间可以互相转换

class Unit:
    def __init__(self,utype:str) -> None:
        self.utype=utype # 单位类型
        self.maps:dict[str,float]={} # 单位与转换系数之间的映射
        self.value:float
    
    @property
    def stand(self)->str:
        return self.units[0]

    @property
    def units(self)->list[str]:
        return list(self.maps.keys())
    
    @property
    def coeff(self)->list[float]:
        return list(self.maps.values())
    
    def set(self,unit:str,value:float):
        assert self.include(unit),f'不符合的单位{unit}'
        index=self.units.index(unit)
        self.value=value/self.coeff[index]
        return self
    
    def get(self,unit:str): # 将指定单位的数值转为标准单位
        assert self.include(unit),f'不符合的单位{unit}'
        index=self.units.index(unit)
        return self.value*self.coeff[index]
    
    def include(self,unit:str):
        if unit in self.units:
            return True
        return False

class ChargeUnit(Unit):
    def __init__(self) -> None:
        super().__init__('charge')
        self.maps={
            'e':1.,
            'c':1.602e-19, # 库伦
        }

class LengthUnit(Unit):
    def __init__(self) -> None:
        super().__init__('length')
        self.maps={
            'bohr':1.,
            'angstrom':0.52918,
            'a':0.52918,
            'nm':0.1
        }
        
class EnergyUnit(Unit):
    def __init__(self) -> None:
        super().__init__('energy')
        self.maps={
            'hartree':1.,
            'ev':27.211,
            'kj/mol':2625,
            'j':4.36e-18
        }

units:list[Unit]=[ChargeUnit(),LengthUnit(),EnergyUnit()]

def trans(unit0:str,unit1:str,value:float):
    """将单位0转换为单位1"""
    unit0=unit0.lower()
    unit1=unit1.lower()
    for unit in units:
        if unit.include(unit0):
            assert unit.include(unit1),f'不符合的单位{unit1}'
            return unit.set(unit0,value).get(unit1)
    raise ValueError(f'不符合的单位{unit0}')

def coeff(unit0:str,unit1:str):
    """获取从单位0到单位1的转换系数"""
    for unit in units:
        if unit.include(unit0):
            assert unit.include(unit1),f'不符合的单位{unit1}'
            return unit.maps[unit1]/unit.maps[unit0]
    raise ValueError(f'不符合的单位{unit0}')