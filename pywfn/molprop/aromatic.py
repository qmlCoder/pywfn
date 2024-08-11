"""
计算分子的芳香性
使用pi键级的标准差表示
"""
from pywfn.base import Mol
from pywfn.atomprop import charge
from pywfn.bondprop import order
import numpy as np

class Calculator:
    def __init__(self, mol:Mol):
        self.mol = mol
        self.caler=order.Calculator(mol)
    
    def calculate(self,ratio:float=0.5):
        result=self.caler.piOrder()
        orders=result[:,-1]
        mean=np.mean(orders)
        stds=np.std(orders)
        return ratio*mean-(1-ratio)*stds