from functools import lru_cache

rawData=[
    [ 1, 'H'],[ 2,'He'],[ 3,'Li'],[ 4,'Be'],[ 5, 'B'],
    [ 6, 'C'],[ 7, 'N'],[ 8, 'O'],[ 9, 'F'],[10,'Ne'],
    [11,'Na'],[12,'Mg'],[13,'Al'],[14,'Si'],[15, 'P'],
    [16, 'S'],[17,'Cl'],[18,'Ar'],[19, 'K'],[20,'Ca'],
    [21,'Sc'],[22,'Ti'],[23, 'V'],[24,'Cr'],[25,'Mn'],
    [26,'Fe'],[27,'Co'],[28,'Ni'],[29,'Cu'],[30,'Zn'],
    [31,'Ga'],[32,'Ge'],[33,'As'],[34,'Se'],[35,'Br'],
    [36,'Kr'],[37,'Rb'],[38,'Sr'],[39, 'Y'],[40,'Zr'],
    [41,'Nb'],[42,'Mo'],[43,'Tc'],[44,'Ru'],[45,'Rh'],
    [46,'Pd'],[47,'Ag'],[48,'Cd'],[49,'In'],[50,'Sn'],
    [51,'Sb'],[52,'Te'],[53, 'I'],[54,'Xe'],[55,'Cs'],
    [56,'Ba'],[57,'La'],[58,'Ce'],[59,'Pr'],[60,'Nd'],
]

class Element:
    def __init__(self,each:list):
        idx,symbol=each
        self.idx=int(idx)
        self.charge=self.idx
        self.symbol=symbol
        self.radius=1.889*0.8

# 因为全局只需要一个Elements实例，因此可以使用类方法
class Elements:
    def __init__(self) -> None:
        self.elements:list[Element]=[Element(each) for each in rawData] #原子列表，一定是按照顺序排列的
        self.symbols:list[str]=[e.symbol for e in self.elements]
        self.charges:list[str]=[e.charge for e in self.elements]
        self.i2e={}
        self.s2e={}

    def __getitem__(self,key:int|str) -> Element:
        if isinstance(key,int):
            assert key in self.charges,f'没有{key}对应元素'
            if key not in self.i2e.keys():
                for element in self.elements:
                    if element.idx==key:
                        self.i2e[key]=element
            return self.i2e[key]
        elif isinstance(key,str):
            assert key in self.symbols,f'没有{key}对应元素'
            if key not in self.s2e.keys():
                for element in self.elements:
                    if element.symbol==key:
                        self.s2e[key]=element
            return self.s2e[key]
        else:
            raise ValueError(f'key必须为str或int类型')

elements=Elements()