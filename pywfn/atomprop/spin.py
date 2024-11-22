"""
此脚本用来计算mulliken自旋
自旋是α电子-β电子
可以通过改变轨道的占据情况来分别计算α和β电子数
"""

from pywfn.base import Mol,Atom
from pywfn.utils import printer
from pywfn.atomprop import charge, lutils
from pywfn.atomprop.charge import Chrgs
from pywfn.maths import CM2PM


from typing import Literal
import numpy as np

class Calculator():
    def __init__(self,mol:"Mol") -> None:
        self.mol=mol
        
    def spin(self,chrg:str='mulliken')->np.ndarray:
        """计算所有原子的自旋"""
        if not self.mol.open: #闭壳层自旋肯定为0
            return np.zeros(self.mol.atoms.natm)
        nmat=self.mol.CM.shape[0] # 系数矩阵行数，基函数数量
        occs_old=self.mol.obtOccs # 记录原本的占据情况
        
        a_occs=occs_old.copy() # 当你需要修改一个变量的时候，
        b_occs=occs_old.copy()
        a_occs[nmat:]=[False]*nmat
        b_occs[:nmat]=[False]*nmat
        elects=[]
        for occs in (a_occs,b_occs):
            self.mol.props.set('obtOccs',occs)
            caler=charge.Calculator(self.mol)
            elect=caler.charge(chrg) # 计算电荷分布
            elects.append(elect)
        a_elect,b_elect=elects
        print(f'atm:{"Na":>10},{"Nb":>10}')
        for i,(ela,elb) in enumerate(zip(a_elect,b_elect)):
            print(f'{i+1:>3d}:{ela:>10.4f},{elb:>10.4f}')
        print('-'*25)
        # 恢复分子属性
        self.mol.props.set('obtOccs',occs_old)
        return -(a_elect-b_elect)
    
    def onShell(self):
        while True:
            printer.options('原子自旋',{
                '1':'Mulliken  电子自旋',
                '2':'Lowdin    电子自旋',
                '3':'Hirshdeld 电子自旋',
            })
            opt=input('请输入自旋类型:')
            match opt:
                case '1':
                    spins=self.spin(chrg='mulliken')
                    for i,spin in enumerate(spins):
                        print(f'{i+1:>3d}:{spin:>10.4f}')
                    print(f'自旋之和: {sum(spins):>10.4f}')
                case '2':
                    spins=self.spin(chrg='lowdin')
                    for i,spin in enumerate(spins):
                        print(f'{i+1:>3d}:{spin:>10.4f}')
                    print(f'自旋之和: {sum(spins):>10.4f}')
                case '3':
                    spins=self.spin(chrg='hirshfeld')
                    for i,spin in enumerate(spins):
                        print(f'{i+1:>3d}:{spin:>10.4f}')
                    print(f'自旋之和: {sum(spins):>10.4f}')
                case _:
                    break