"""
定义分子结构编辑器
可以用来生成刚性扫描的结构
"""

from pywfn.base import Mole
from pywfn.maths import points_rotate

import numpy as np

class Editor:
    def __init__(self,mol:Mole) -> None:
        self.mol=mol
    
    def search_group(self,start:int,excludes:list[int]):
        """
        搜索原子组，从指定原子开始用深度搜索算法查找原子组
        start:开始的原子
        eccludes:排除的原子
        """
        group=[start]
        natm=self.mol.atoms.num
        for i in range(natm):
            if i==len(group):break
            atm=group[i]
            nebs=self.mol.atom(atm).neighbors
            # print(nebs)
            for neb in nebs:
                if neb in group:continue
                if neb in excludes:continue
                group.append(neb)
        return group
    
    def rotate_bond(self,atm1:int,atm2:int,angle:float):
        """
        旋转两个原子之间的键
        """
        # 搜索原子组
        group=self.search_group(atm1,excludes=[atm2])
        coords=self.mol.coords.copy()
        gidxs=[e-1 for e in group]
        points=coords[gidxs,:]
        center=self.mol.atom(atm1).coord
        axis=self.mol.atom(atm2).coord-center
        result=points_rotate(points,center,axis,angle)
        coords[gidxs,:]=result
        syms=self.mol.geome.syms.copy()
        self.mol.geome.build(syms,coords)
        return self.mol
    
    def add_ringBq(self): # 在环中心添加Bq原子
        pass
        