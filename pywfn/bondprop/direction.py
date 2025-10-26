"""
计算与键相关的方向
"""
from pywfn.base.mole import Mole
from pywfn.maths import points_rotate
import numpy as np

class Calculator:
    def __init__(self,mole:Mole) -> None:
        self.mole=mole

    def verticals(self,atm1:int,atm2:int)->np.ndarray:
        """计算垂直于键轴一圈的方向

        Args:
            atm1 (int): 原子1索引
            atm2 (int): 原子2索引

        Returns:
            np.ndarray: 方向数组[35,3]
        """
        bond=self.mole.bonds.get(atm1,atm2)
        assert bond is not None,"No bond between atoms"
        bdir=bond.vector # 键轴向量
        anyv=np.random.rand(3) # 随机生成向量

        bdir/=np.linalg.norm(bdir)
        anyv/=np.linalg.norm(anyv)

        vect=np.cross(bdir,anyv) # 垂直于键轴的向量
        cent=self.mole.atom(atm1).coord
        point=cent+vect
        point=point.reshape(1,3)

        angs=np.linspace(0,np.pi*2,35,endpoint=False)
        result=np.zeros(shape=(len(angs),3))
        for a,ang in enumerate(angs):
            result[a]=points_rotate(point,cent,bdir,ang)-cent
        
        return result
