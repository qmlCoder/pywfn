"""
该脚本用于测试给pi键级加方向矫正
"""

import numpy as np
from pywfn.atomprop import direction
from pywfn.base import Mole
from pywfn.bondprop import order


def get_ang(dir1, dir2) -> float:
    # 1. 计算向量点积
    dir1 = dir1 / np.linalg.norm(dir1)
    dir2 = dir2 / np.linalg.norm(dir2)
    dot_product = np.dot(dir1, dir2)
    # 2. 计算向量的模长
    norm_a = np.linalg.norm(dir1)
    norm_b = np.linalg.norm(dir2)
    assert norm_a != 0 and norm_b != 0
    # 3. 计算余弦值（限制范围[-1,1]，避免浮点误差报错）
    cos_theta = np.clip(dot_product / (norm_a * norm_b), -1.0, 1.0)
    # 4. 计算夹角（弧度）
    ang = np.arccos(cos_theta)
    # print(dir1, dir2, ang)
    return ang


path = r"E:\code\articles\MOCV\moles\B12H12_21_wfn.log"
mole = Mole.from_file(path)


dir_caler = direction.Calculator(mole)
natm = mole.atoms.len()

stms = {}
for atm in range(natm):
    stm = dir_caler.LCS(atm)
    stms[atm] = stm


caler = order.Calculator(mole)
stms, atos, omat = caler.pi_mocv()


# 两个B-H键向量之间的夹角
def get_fmat1():
    fmat = np.zeros((natm, natm))
    for i in range(natm):
        for j in range(natm):
            if not mole.atoms.is_bond(i, j):
                continue
            if i == j:
                continue
            if i not in stms or j not in stms:
                continue
            diri = np.array(stms[i].ex)
            dirj = np.array(stms[j].ex)
            ang = get_ang(diri, dirj)
            val = np.cos(ang)
            fmat[i, j] = val**0.5
    return fmat


# B-H键与B-B键法向量之间夹角的两个cos值的乘积
def get_fmat2():
    fmat = np.zeros((natm, natm))
    for i in range(natm):
        for j in range(natm):
            if not mole.atoms.is_bond(i, j):
                continue
            if i == j:
                continue
            if i not in stms or j not in stms:
                continue
            stma = dir_caler.LCS(i, j)
            stmb = dir_caler.LCS(j, i)
            dira = np.array(stma.ex)
            dirb = np.array(stmb.ex)
            diri = np.array(stms[i].ex)
            dirj = np.array(stms[j].ex)
            angi = get_ang(dira, diri)
            angj = get_ang(dirb, dirj)
            print(i + 1, j + 1, angi / np.pi * 180)
            val = np.cos(angi) * np.cos(angj)
            fmat[i, j] = val**0.5
    return fmat


# 两个B-H键向量之间的夹角
def get_fmat3():
    fmat = np.zeros((natm, natm))
    for i in range(natm):
        for j in range(natm):
            if not mole.atoms.is_bond(i, j):
                continue
            if i == j:
                continue
            if i not in stms or j not in stms:
                continue
            diri = np.array(stms[i].ex)
            dirj = np.array(stms[j].ex)
            ang = get_ang(diri, dirj)
            val = np.cos(np.pi / 2 - ang)
            fmat[i, j] = val**0.5
    return fmat


# fmat = get_fmat1()  # 0.1433
fmat = get_fmat2()  # 0.1822
# fmat = get_fmat3()  # 0.2027

omat *= fmat

for i in range(natm):
    for j in range(natm):
        if not mole.atoms.is_bond(i, j):
            continue
        if i == j:
            continue
        val = omat[i, j]
        print(f"omat[{i:>2}, {j:>2}] = {val}")
