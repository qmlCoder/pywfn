# 用来管理程序的一些设置
from pathlib import Path
import numpy as np

BOND_LIMIT=1.7 # 判断量原子之间是否成键的长度限制

IF_DEBUG=True # 是否开启debug,控制项目中所有的打印,避免与shell中的print冲突
IF_SHELL=False # 是否在shell中执行

BASE_VECTOR=np.array([0,0,1]) # 标准向量，求原子法向量和轨道方向的时候要与该向量夹角小于90°

IF_ORBITAL_ORDER=False #是否计算键级中每个轨道的成分，尤其是mayer键级拆分成每个轨道的成分很麻烦

RENDER_ATOM_RANGE=20  # 原子轨道的渲染范围

RENDER_CLOUD_STEP=0.25  # 渲染间隔

RENDER_CLOUD_BORDER=3

BOHR_RADIUS=1.889

TEMPLATES:dict[str,Path]={
    'gif':None,
    'si':None
}