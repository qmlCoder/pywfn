import sys
sys.path.append("D:\code\pywfn")

from pywfn.base import Mol

# 读取mol文件
from pywfn.reader import MolReader
path="D:\R-C4-OH-.mol"
reader=MolReader(path)
result=reader.read_coords()
mol=Mol(reader=reader)
print(mol.coords)