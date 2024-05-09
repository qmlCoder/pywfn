"""
定义一些常见的自定义类型，否则在注释的时候都无法区分
"""
import numpy as np

class Coord(np.ndarray):
    def __init__(self,xyz:list[float]) -> None:
        super().__init__(xyz)
        self.dtype=float

class Lmn(np.ndarray):
    def __init__(self,lmn:list[float]) -> None:
        super().__init__(lmn)
        self.dtype=int