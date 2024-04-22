import sys
sys.path.append("D:\code\pywfn")

from pywfn.base import Mol
from pywfn.reader import LogReader
from pywfn.writer import cubWriter

path="D:\BaiduSyncdisk\Articles\HFV\gfile\CH4.log"
reader=LogReader(path)

mol=Mol(reader=reader)

writer=cubWriter(mol)
writer.obts=[0,1,2]
writer.save('CH4_0,1,2')
print(writer.filePath)