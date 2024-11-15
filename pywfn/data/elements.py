"""
记录元素数据
- 元素序号
- 元素符号
- 原子半径 https://chemistry-europe.onlinelibrary.wiley.com/doi/full/10.1002/chem.200800987
"""

from functools import lru_cache

rawData = [
    [ 0, "Bq", 25],
    [ 1, "H" , 32],
    [ 2, "He", 46],
    [ 3, "Li",133],
    [ 4, "Be",102],
    [ 5, "B" , 85],
    [ 6, "C" , 75],
    [ 7, "N" , 71],
    [ 8, "O" , 63],
    [ 9, "F" , 64],
    [10, "Ne", 67],
    [11, "Na",155],
    [12, "Mg",139],
    [13, "Al",126],
    [14, "Si",116],
    [15, "P" ,111],
    [16, "S" ,103],
    [17, "Cl", 99],
    [18, "Ar", 96],
    [19, "K" ,196],
    [20, "Ca",171],
    [21, "Sc",148],
    [22, "Ti",136],
    [23, "V" ,134],
    [24, "Cr",122],
    [25, "Mn",119],
    [26, "Fe",116],
    [27, "Co",111],
    [28, "Ni",110],
    [29, "Cu",112],
    [30, "Zn",118],
    [31, "Ga",124],
    [32, "Ge",121],
    [33, "As",121],
    [34, "Se",116],
    [35, "Br",114],
    [36, "Kr",117],
    [37, "Rb",210],
    [38, "Sr",185],
    [39, "Y" ,163],
    [40, "Zr",154],
    [41, "Nb",147],
    [42, "Mo",138],
    [43, "Tc",128],
    [44, "Ru",125],
    [45, "Rh",125],
    [46, "Pd",120],
    [47, "Ag",128],
    [48, "Cd",136],
    [49, "In",142],
    [50, "Sn",140],
    [51, "Sb",140],
    [52, "Te",136],
    [53, "I" ,133],
    [54, "Xe",131],
    [55, "Cs",232],
    [56, "Ba",196],
    [57, "La",180],
    [58, "Ce",163],
    [59, "Pr",176],
    [60, "Nd",174],
]


class Element:
    def __init__(self, idx:int,sym:str,rad:float):
        self.idx = idx
        self.sym = sym
        self.rad = rad/100/0.529177
        self.symbol = self.sym
        self.charge = self.idx
        self.atomic = self.idx
        self.radius = self.rad


# 因为全局只需要一个Elements实例，因此可以使用类方法
class Elements:
    def __init__(self) -> None:
        self.elements: list[Element] = [
            Element(idx,sym,rad) for idx,sym,rad in rawData
        ]  # 原子列表，一定是按照顺序排列的
        self.symbols: list[str] = [e.symbol for e in self.elements]
        self.charges: list[int] = [e.charge for e in self.elements]
        self.i2e = {}
        self.s2e = {}

    def __getitem__(self, key: int | str) -> Element:
        if isinstance(key, int):
            assert key in self.charges, f"没有{key}对应元素"
            if key not in self.i2e.keys():
                for element in self.elements:
                    if element.idx == key:
                        self.i2e[key] = element
            return self.i2e[key]
        elif isinstance(key, str):
            assert key in self.symbols, f"没有{key}对应元素"
            if key not in self.s2e.keys():
                for element in self.elements:
                    if element.symbol == key:
                        self.s2e[key] = element
            return self.s2e[key]
        else:
            raise ValueError(f"key必须为str或int类型")


elements = Elements()
