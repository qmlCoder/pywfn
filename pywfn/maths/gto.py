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

分子轨道波函数
原子轨道波函数
原子波函数->原子电子密度
"""

import numpy as np
from numpy import sqrt, ndarray

from pywfn import base
from pywfn.utils import printer
from pywfn.data import Basis
from pywfn.maths import flib

π = np.pi
e = np.e


def fac2(num):  # 该例中一定是奇数
    """计算双阶乘"""
    # num+=1
    if num <= 1:
        fac = 1
    else:
        fac = np.prod(np.arange(1, num + 1, 2))
    return fac


class Gto:
    def __init__(self, mol: "base.Mol"):
        self.mol = mol
        self.basis = mol.basis
        self.angs = [0, 1, 2, 3, 4]  # 默认全部绘制

    # def ato(self, grid: ndarray, atm: int, obts: list[int]):
    #     """
    #     计算原子的高斯型波函数数值
    #     grid:格点，以原子为中心
    #     atm:原子
    #     obt:轨道
    #     """
    #     assert len(grid.shape) == 2, f"pos的维度应该为2"
    #     assert grid.shape[1] == 3, f"pos的形状应为[n,3]，当前为{grid.shape}"
    #     R2 = np.sum(grid**2, axis=1)  # x^2+y^2+z^2

    #     atom = self.mol.atom(atm)
    #     u, l = atom.obtBorder

    #     atms = self.mol.obtAtms[u:l]
    #     shls = self.mol.obtShls[u:l]
    #     syms = self.mol.obtSyms[u:l]
    #     lmns = [Basis.sym2lmn(sym) for sym in syms]
    #     coef = self.mol.CM[u:l, obts]  # 二维矩阵
    #     expl = []
    #     coel = []
    #     ncsl = []  # 记录每个收缩轨道的大小
    #     for i in range(len(atms)):  # 该原子的行索引
    #         lmn = lmns[i]
    #         atm = atms[i]
    #         shl = shls[i]
    #         ang = sum(lmn)
    #         basis = self.mol.basis.get(atom.atomic, shl, ang)
    #         exps = [b.exp for b in basis]
    #         coes = [b.coe for b in basis]
    #         expl += exps
    #         coel += coes
    #         ncsl.append(len(exps))
    #     wfns = self.agf(expl, coel, ncsl, lmns, R2, grid, coef, obts)
    #     return wfns

    # def agf(
    #     self,
    #     expl: list[float],
    #     coel: list[float],
    #     ncsl: list[int],
    #     lmns: list[int],
    #     R2: np.ndarray,
    #     pos: np.ndarray,
    #     coef: np.ndarray,
    #     obts: list[int],
    # ) -> np.ndarray:
    #     nexp = len(expl)
    #     j = 0
        
    #     wfns = np.zeros(len(pos), dtype=np.float32)
    #     for i in range(len(lmns)):
    #         nc = ncsl[i]  # 收缩大小
    #         lmn = lmns[i]
    #         coes = coel[j : j + nc]
    #         exps = expl[j : j + nc]
    #         j += nc
    #         wfn = self.cgf(exps, coes, lmn, R2, pos) # 每一个原子轨道都是提前定义好的不变的
    #         for obt in obts:
    #             wfns += coef[i, obt] * wfn
    #     return wfns
    
    def cgf(self,exps: list[float],coes: list[float],lmn: list[int],grids: np.ndarray,coord:np.ndarray) -> np.ndarray:
        # wfn = self.cgf_prot(exps,coes,lmn,grid)
        wfn = self.cgf_fort(exps,coes,lmn,grids,coord)
        return wfn

    def cgf_prot(self,exps: list[float],coes: list[float],lmn: list[int],grid: np.ndarray,coord:np.ndarray) -> np.ndarray:
        """
        收缩高斯函数，线性组合之前的波函数
        pos:坐标[n,3]，以原子为中心
        c1*f1+c2*f2...
        """
        wfns = np.zeros(len(grid))
        for exp, coe in zip(exps, coes):
            wfns += coe * self.gtf(exp, grid, lmn, coord)
        return wfns
    
    def cgf_fort(self,exps: list[float],coes: list[float],lmn: list[int],grids: np.ndarray,coord:np.ndarray) -> np.ndarray:
        """收缩基函数"""
        cmax=len(exps)
        nc=len(exps)
        alps=np.array(exps)
        coes=np.array(coes)
        ngrid=grids.shape[0]
        l,m,n=lmn
        wfns=flib.cgf(cmax,nc,alps,coes,ngrid,grids,coord,l,m,n)
        return wfns

    def gtf(self, exp: float, grids: ndarray, lmn: list[int], coord:np.ndarray) -> ndarray:
        """
        计算指定点gto函数的值
        exp:高斯指数
        grids:格点[n,3]，以原子为中心
        lmn:角动量决定的x,y,z指数
        coord:原子坐标
        """
        # i=0,(2i-1)=-1,(2i-1)!!=1
        # i=1,(2i-1)= 1,(2i-1)!!=1
        # i=2,(2i-1)= 3,(2i-1)!!=3
        # print('gtf',pos[0,:])
        facs=[1,1,3] # n与阶乘的对应关系？
        l, m, n = lmn
        x, y, z = (grids-coord).T # 以原子坐标为中心
        R2 = np.sum(grids**2, axis=1)  # x^2+y^2+z^2
        ang = sum(lmn)  # 角动量
        fac = facs[l] * facs[m] * facs[n]  # 双阶乘
        N = (2 * exp / π) ** (3 / 4) * np.sqrt((4 * exp) ** ang / fac)
        wfn = x**l * y**m * z**n * np.exp(-exp * R2) * N
        return wfn
