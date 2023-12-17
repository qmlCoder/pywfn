"""
轨道挑选法+HMO键级公式
"""
import numpy as np

from pywfn.base import Mol
from pywfn.bondorder import lutils,Caler

class Calculator(Caler):
    def __init__(self,mol:Mol) -> None:
        self.mol=mol
        self.oe=self.mol.oE
        self.bond:list[int]=None

    def calculate(self):
        idx1,idx2=self.bond
        cenAtom=self.mol.atom(idx1)
        aroAtom=self.mol.atom(idx2)
        obts=self.mol.O_obts

        As=self.mol.As[obts]
        piUnits=[lutils.judgeOrbital(cenAtom,aroAtom,obt,cenAtom.get_Normal()) for obt in obts]
        cenRes=np.sqrt(np.sum(cenAtom.OC[:,obts]**2,axis=0))/As
        aroRes=np.sqrt(np.sum(aroAtom.OC[:,obts]**2,axis=0))/As
    
        orders=cenRes*aroRes*np.array(piUnits)*self.mol.oE
        lutils.printOrders(orders,[obt+1 for obt in obts])
        return sum(orders)