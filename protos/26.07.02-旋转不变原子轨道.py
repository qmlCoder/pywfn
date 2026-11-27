"""
找到随x轴旋转保持不变的局部坐标系下的原子轨道组合
"""

import functools
import itertools
import math

import numpy as np
from pywfn.datas import mocv


def get_lmat():
    ex = np.array([1, 0, 0])
    ez = np.array([0, np.random.rand(), np.random.rand()])
    ez = ez / np.linalg.norm(ez)
    ey = np.cross(ez, ex)
    ey = ey / np.linalg.norm(ey)
    return np.array([ex, ey, ez]).T


print(get_lmat())
print(get_lmat())

# rcoefs = np.random.rand(6, 1)  # 原始系数
# for i in range(10):
#     lmat = get_lmat()
#     tmat = mocv.tmat_d(lmat)
#     lcoefs = tmat @ rcoefs  # 局部坐标系下的系数
#     lcoefs = lcoefs * np.array([0, 0, 0, 1, 1, 0]).reshape(6, 1)
#     fcoefs = np.linalg.inv(tmat) @ lcoefs
#     print(fcoefs.flatten())


def find_didx():
    rcoefs = np.random.rand(6, 1)  # 原始系数
    fix = np.zeros((6, 1))
    for idxs in itertools.combinations(range(6), 2):
        fix[idxs, 0] = 1
        data = np.zeros((5, 6))
        for i in range(5):
            lmat = get_lmat()
            tmat = mocv.tmat_d(lmat)
            lcoefs = tmat @ rcoefs  # 局部坐标系下的系数
            lcoefs = lcoefs * fix
            fcoefs = np.linalg.inv(tmat) @ lcoefs
            data[i] = fcoefs.flatten()
        find = np.allclose(data, data[0, :])
        if find:
            print(idxs)


def find_fidx():
    rcoefs = np.random.rand(10, 1)  # 原始系数
    fix = np.zeros((10, 1))
    allows = [4, 5, 3, 9]
    for idxs in itertools.combinations(allows, 5):
        fix[idxs, 0] = 1
        data = np.zeros((5, 10))
        for i in range(5):
            lmat = get_lmat()
            tmat = mocv.tmat_f(lmat)
            Mi = tmat.T
            Mt = np.linalg.inv(Mi)
            lcoefs = Mt @ rcoefs  # 局部坐标系下的系数
            lcoefs = lcoefs * fix
            fcoefs = Mi @ lcoefs  # 全局坐标系下的系数
            data[i] = fcoefs.flatten()
        find = np.allclose(data, data[0, :])
        if find:
            print(data)
            print(idxs)


find_didx()
print("-" * 10)
find_fidx()
