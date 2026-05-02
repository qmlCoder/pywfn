"""
计算原子的反应活性，福井函数、parr函数等，一般需要多个分子才行
"""
from pywfn.base import Mole
from pywfn.atomprop import charge
import numpy as np
from pywfn import config
from pywfn import _core


class Calculator:
    def __init__(self, mole: Mole) -> None:
        self.mole = mole
        self.core = _core.atomprop.activity.Calculator( # type: ignore
            mole.core)  # type: ignore # 核心计算器

    # 福井函数
    def fukui(self, mole_n: Mole, mole_p: Mole, ctype: str = "mulliken"):
        """计算所有原子的福井函数

        Args:
            molN (Mole): 负电荷分子
            molP (Mole): 正电荷分子
            ctype (Chrgs, optional): 电荷类型. Defaults to 'mulliken'.

        Returns:
            np.ndarray: 福井函数[n,7](N,N+1,N-1,f-,f+,f0,df)
        """
        vals = self.core.fukui(mole_n.core, mole_p.core, ctype)
        return np.array(vals)

    def fukui_pi(self, mole_n: Mole, mole_p: Mole, ctype: str = 'mulliken'):
        """基于pi电子数的福井函数"""
        dirs, vals = self.core.fukui_pi(mole_n.core, mole_p.core, ctype)
        return dirs, np.array(vals)
    
    # 某个原子的方向福井函数
    def fukui_dir(self, atm: int, dir: list[float], mole_n: Mole, mole_p: Mole, ctype: str = "mulliken"):
        """计算方向福井函数

        Args:
            atms (int): 要计算的原子
            dirs (np.ndarray): 原子的方向，可以为多个[n,3]
            molN (Mole): 多一个电子的分子
            molP (Mole): 少一个电子的分子
            ctype(str): 电荷类型,默认为mulliken

        Returns:
            np.ndarray: 指定方向的福井函数[n,6](N,N+1,N-1,f-,f+,f0,df)
        """
        vals = self.core.fukui_dir(
            atm, dir, mole_n.core, mole_p.core, ctype)
        return np.array(vals)

    # 自由价
    def freev(self, otype:str="mayer"):
        """计算指定原子的自由价

        Args:
            atm (int): 原子索引

        Returns:
            np.ndarray: 自由价[d](val)
        """
        return self.core.freev(otype)

    def freev_pi(self,otype:str="mayer"):
        """使用pi键级计算的自由价

        Args:
            otype (str, optional): 键级类型 order type. Defaults to "mayer".
        """
        dirs,vals=self.core.freev_pi(otype)
        return dirs,np.array(vals)

    def freev_dir(self,atm:int,dir:list[float],nebs:list[int],otype:str="mayer"):
        return self.core.freev_dir(atm,dir,nebs,otype)


    
