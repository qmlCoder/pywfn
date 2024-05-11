"""
实现单位转换
只需要实现化学中常用的单位转换就行
根据单位的类型，将单位转换为标准值，再将标准值转换为需要的单位
可以直接转换数值，也可以获得转换系数
"""

# 定义支持的单位和类型，同一个类型的单位之间可以互相转换
units={
    'energy':[''],
    'length':['angstrom','borh']
}

def trans(unit0:str,unit1:str,value:float):
    """将单位0转换为单位1"""
    utype=''
    if not utype:return None