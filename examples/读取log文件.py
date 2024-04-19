import sys
sys.path.append("D:\code\pywfn")

from pywfn.base import Mol
from pywfn.reader import LogReader

path="D:\BaiduSyncdisk\Articles\HFV\gfile\CH4.log"
print(path[-4:])
reader=LogReader(path)

occs=reader.get_obtOccs()
print(occs)