"""
休克尔分子轨道法轨道的可视化

将空间格点转换为原子局部坐标系的点，计算这些点上的波函数值
"""

from pywfn.base.mole import Mole


class Calculator:
    def __init__(self, mole: Mole):
        self.mole = mole
