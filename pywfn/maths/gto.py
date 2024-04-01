"""
传入一堆位置点，计算在这些位置的各种函数数值
传入的数组大小为[n,3],一行分别是x,y,z
基组的种类是有限的吗？有没有通用的表达形式？为每一个基组都定义不同的计算类型
将多次的向量计算转为多次的矩阵运算
坐标数量 N
收缩数量 S
原子轨道数量 L
原子数量 N
一个gto可以分为两部分

一个原子轨道由许多电子轨道组成，每个电子轨道最多容纳两个电子
每个电子轨道由多个gto组成，用多个gto去拟合(线性组合)一个电子轨道

每一个电子轨道有一个对应的分子轨道系数

电子轨道又可以分为不同层

每一层有不同的角动量 l+m+n
每一个角动量又有三个分量 l,m,n
每个电子轨道对应一个角动量分量
组成电子轨道的gto的系数相同但指数不同c1*gto(a)+c2*gto(a)+c3*gto(a)

代码的可读性和速度如何取舍？ 我还是选择了后者

可以将gto的计算放到多个进程中进行
"""
import numpy as np
from numpy import sqrt,ndarray

from pywfn import base
from pywfn.utils import printer

π=np.pi
e=np.e

class Gto:
    def __init__(self,mol:"base.Mol"):
        self.mol=mol
        self.basis=mol.basis
        self.angs=[0,1,2,3,4] # 默认全部绘制

    def bind(self,atom:int,obt:int):
        atom:"base.Atom"=self.mol.atom(atom)
        datas=self.basis.get(atom.atomic)
        coeff=atom.obtCoeffs[:,obt]
        data_=[]
        idxo=0
        coes={}
        for shl,ang,exp,coe in datas:
            for l,m,n in self.basis.lmn(ang):
                key=f'{shl+1}{self.basis.lName(l,m,n)}'
                if key not in coes.keys(): # 保证唯一性
                    coes[key]=coeff[idxo]
                    idxo+=1
                data_.append([coes[key],exp,coe,l,m,n])
        assert idxo==len(coeff),f'系数没有匹配完全{idxo=},{coeff=}'
        return data_
    
    def agto(self,pos:ndarray,atom:int,obt:int):
        assert pos.shape[1]==3,f'pos的形状应为[n,3]，当前为{pos.shape}'
        pos*=1.889 #这一步很关键，将埃转为波尔
        R=np.sum(pos**2,axis=1)
        values=np.zeros(len(pos),dtype=np.float32)
        data_=self.bind(atom,obt)
        for C,exp,coe,l,m,n in printer.track(data_,f'agto:atom={atom}'):
            if not (l+m+n in self.angs):continue #可以筛选角动量
            if abs(C)<1e-6:continue #系数很小的忽略
            values+=C*self.gto(exp,coe,R,pos,l,m,n)
        return values

    def gto(self,exp:float,coe:float,R:ndarray,pos:ndarray,l:int,m:int,n:int)->ndarray:
        """
        计算指定点gto函数的值
        c:收缩系数
        a:高斯指数
        pos:坐标
        R:坐标到原点距离平方
        lmn:角动量决定的x,y,z指数
        """
        x,y,z=pos.T
        v:ndarray=coe*(2*exp/π)**(3/4)*2*sqrt(exp)* x**l * y**m * z**n * e**(-exp*R)
        return v