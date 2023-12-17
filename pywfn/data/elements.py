from functools import lru_cache

rawData=[
    [1,'H'],[2,'He'],[3,'Li'],[4,'Be'],[5,'B'],
    [6,'C'],[7,'N'],[8,'O'],[9,'F'],[10,'Ne'],
    [11,'Na'],[12,'Mg'],[13,'Al'],[14,'Si'],[15,'P'],
    [16,'S'],[17,'Cl'],[18,'Ar'],[19,'K'],[20,'Ca'],
]

class Element:
    def __init__(self,each:list):
        idx,symbol=each
        self.idx=int(idx)
        self.charge=self.idx
        self.symbol=symbol
        self.radius=0.8

class Elements:
    def __init__(self) -> None:
        self.elements:list[Element]=[Element(each) for each in rawData] #原子列表，一定是按照顺序排列的

    @lru_cache
    def __getitem__(self,key:int|str) -> Element:
        if isinstance(key,int):
            for element in self.elements:
                if element.idx==key:
                    return element
            raise ValueError(f'没有{key}对应元素')
        elif isinstance(key,str):
            for element in self.elements:
                if element.symbol==key:
                    return element
            raise ValueError(f'没有{key}对应元素')
        else:
            raise ValueError(f'key必须为str或int类型')

elements=Elements()