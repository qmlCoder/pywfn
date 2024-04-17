"""
计算波函数值
"""
import sys
sys.path.append('D:\code\pywfn')
from pywfn.base import Mol
from pywfn.reader import LogReader
import numpy as np

path="D:\BaiduSyncdisk\gfile\elements\H2.out"
mol=Mol(reader=LogReader(path))
mol.bohr=True
pos=np.array([
    [-3.936929,-3.936929,-4.566838]
])
val=mol.get_wfnv(1,pos,[1,2])
print(val)