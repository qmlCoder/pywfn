import matplotlib.pyplot as plt
import numpy as np
from pywfn.base import Mole

root = r"E:\code\articles\Research Diary\B12H12 Study\moles"

mole = Mole.from_file(rf"{root}\B4H4.log")
frag1 = Mole.from_file(rf"{root}\B4H4-f1.log")
frag2 = Mole.from_file(rf"{root}\B4H4-f2.log")

A = mole.coefs.get_cmat()
B1 = frag1.coefs.get_cmat()
B2 = frag2.coefs.get_cmat()

print(A.shape, B1.shape, B2.shape)

nato = A.shape[0]
nato1 = B1.shape[0]
nato2 = B2.shape[0]

B = np.zeros((nato, nato))

B[:nato1, :nato1] = B1
B[nato1:, nato1:] = B2


C = np.linalg.inv(B) @ A

S = mole.basis.get_smat()

S_fmo = B.T @ S @ B  #  FMO 之间的重叠积分

# print(C.T @ S_fmo @ C) # 验证 C @ S_fmo @ C 的结果

W = np.zeros_like(C)

# for i in range(C.shape[1]):
#     ci = C[:, i]  # 每个分子轨道中，片段轨道的系数
#     W[:, i] = ci * (S_fmo @ ci)  # 每一个分子轨道中，片段轨道的权重
# print(W.sum(axis=0))

C2 = C**2
W = C2 / C2.sum(axis=0)
print(W.sum(axis=0))

print("C")
for i in range(nato):
    if i < nato1:
        print(f"frag1 {i:>3}", end="")
    else:
        print(f"frag2 {(i - nato1):>3}", end="")
    for j in range(nato):
        print(f"{C[i, j]:>10.4f}", end="")
    print()

print("W")
for i in range(nato):
    if i < nato1:
        print(f"frag1 {i:>3}", end="")
    else:
        print(f"frag2 {(i - nato1):>3}", end="")
    for j in range(nato):
        print(f"{W[i, j]:>10.4f}", end="")
    print()
