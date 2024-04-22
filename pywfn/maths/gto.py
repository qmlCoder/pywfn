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

def fac2(num): # 该例中一定是奇数
    """计算双阶乘"""
    if num<=1:
        fac = 1
    else:
        fac = np.prod(np.arange(1,num+1,2))
    return fac

class Gto:
    def __init__(self,mol:"base.Mol"):
        self.mol=mol
        self.basis=mol.basis
        self.angs=[0,1,2,3,4] # 默认全部绘制

    def bind(self,atm:int,obt:int): # 将基组信息与轨道系数绑定在一起
        atom=self.mol.atom(atm)
        basis=self.basis.get(atom.atomic)
        coeff=atom.obtCoeffs[:,obt]
        params=[]
        idxo=0
        coes={}
        for _,shl,ang,exp,coe in basis:
            for l,m,n in self.basis.lmn(ang): # 角动量转为l,m,n
                key=f'{shl}{self.basis.lName(l,m,n)}'
                if key not in coes.keys(): # 保证唯一性
                    coes[key]=coeff[idxo]
                    idxo+=1
                params.append([coes[key],exp,coe,l,m,n]) # 一个角动量可以对应多个轨道系数
        assert idxo==len(coeff),f'系数没有匹配完全{idxo=},{coeff=}'
        return params
    
    def agto(self,pos:ndarray,atm:int,obt:int):
        """
        计算原子的高斯型波函数数值
        pos:坐标
        atm:原子
        obt:轨道
        """
        assert len(pos.shape)==2,f'pos的维度应该为2'
        assert pos.shape[1]==3,f'pos的形状应为[n,3]，当前为{pos.shape}'
        # pos*=1.889 #这一步很关键，将埃转为波尔 !取分子坐标时可转为波尔
        R2=np.sum(pos**2,axis=1) # x^2+y^2+z^2
        values=np.zeros(len(pos),dtype=np.float32)
        params=self.bind(atm,obt)
        for C,exp,coe,l,m,n in params:
            if l+m+n not in self.angs:continue #可以筛选角动量
            if abs(C)<1e-6:continue #系数很小的忽略
            wfn=C*self.gto(exp,coe,R2,pos,l,m,n)
            # assert np.max(np.abs(wfn))<1,'波函数值不合理'
            values+=wfn
        return values

    def gto(self,exp:float,coe:float,R2:ndarray,pos:ndarray,l:int,m:int,n:int)->ndarray:
        """
        计算指定点gto函数的值
        coe:收缩系数
        exp:高斯指数
        pos:坐标[n,3]
        R:坐标到原点距离平方
        lmn:角动量决定的x,y,z指数
        """

        x,y,z=pos.T
        # N=coe*(2*exp/π)**(3/4)*2*sqrt(exp) #归一化系数
        ang=l+m+n # 角动量
        fac=fac2(2*l-1)*fac2(2*m-1)*fac2(2*n-1) # 双阶乘
        N=((4*exp)**ang/fac)**(1/2)*(2/π)**(3/4) # 归一化系数
        wfn= x**l * y**m * z**n * e**(-exp*R2) *N *coe
        return wfn