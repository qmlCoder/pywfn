from pywfn import core
import numpy as np

class Calculator:
    def __init__(self,mole) -> None:
        self.mole=mole
        self.caler=core.moleprop.orbital.Calculator(mole.mole) # type: ignore # 核心计算器


class Deco:
    def __init__(self,base:np.ndarray,skeep:list[int],pkeep:list[int],dkeep:list[int]) -> None:
        self.deco=core.moleprop.orbital.Deco(base,skeep,pkeep,dkeep) # type: ignore # 核心计算器