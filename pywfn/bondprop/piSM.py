"""
轨道挑选法+mayer键级公式
mayer键级需要重叠矩阵,可一次性计算所有键级(此时原子的法向量由三个原子决定,而不是垂直于键轴)
"""
import numpy as np

from pywfn.base import Mol
from pywfn.bondprop import lutils,Caler
from pywfn import maths

class Calculator(Caler):
    def __init__(self,mol:Mol) -> None:
        self.mol=mol
        self.bond:list[int]=None

    def calculate(self)->float:
        """指定两个原子,计算π键键级"""
        idx1,idx2=self.bond
        atom1=self.mol.atom(idx1)
        atom2=self.mol.atom(idx2)
        normal=atom1.get_Normal()
        if normal is None:
            idxn=maths.search_sp2(idx1,self.mol)
            if idxn is None:
                normal=atom1.get_vertObt(atom2.idx,atom1.idx)
            else:
                normal=self.mol.atom(idxn).get_Normal(idx2)

        obts=self.mol.O_obts
        CM_=np.zeros_like(self.mol.CM) # 拷贝一份，然后将不是π轨道的那些变成0

        piObts=[]
        for atom in [atom1,atom2]: # 修改每个原子对应的系数矩阵
            if atom.symbol=='H':continue
            a_1,a_2=atom.obtRange
            for orbital in obts:
                judgeRes=lutils.judgeOrbital(atom1,atom2,orbital,normal)
                if judgeRes==0:continue # 如果是π轨道
                CM_[a_1:a_2,orbital]=self.mol.CM[a_1:a_2,orbital]
                piObts.append(orbital)
        oe=1 if self.mol.isOpenShell else 2
        PM_=lutils.CM2PM(CM_,obts,oe)
        SM=self.mol.SM
        PS=PM_@SM

        a1,a2=atom1.obtRange
        b1,b2=atom2.obtRange
        order=np.sum(PS[a1:a2,b1:b2]*PS[b1:b2,a1:a2])
        piObts=set(piObts)
        piObts=[self.mol.obtStr[o] for o in piObts]
        lutils.formPrint(contents=[piObts],eachLength=10)
        return order