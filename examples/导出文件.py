import sys
sys.path.append("D:\code\pywfn\src")

from pywfn.base import Mol
from pywfn.reader import LogReader,MolReader



# 读入文件
path="D:\R-C4-OH-.mol"
reader=MolReader(path)
mol=Mol(reader=reader)


# from pywfn.writer import cubWriter
# writer=cubWriter(mol)
# writer.obts=[0,1,2]
# writer.save('CH4_0,1,2')
# print(writer.filePath)


# from pywfn.writer import MolWriter
# writer=MolWriter(mol)
# writer.build()
# print(writer.temp)
# print(mol.formula)
# writer.save('CH4.mol')



from pywfn.writer import gjfWriter
writer=gjfWriter(mol)

writer.save('CH4.gjf')
