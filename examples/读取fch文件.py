import sys
sys.path.append("D:\code\pywfn")

from pywfn.base import Mol
from pywfn.reader import FchReader

path="D:\BaiduSyncdisk\Articles\HFV\gfile\CH4.fch"
reader=FchReader(path)

reader.search_title()
reader.read_atoms()
reader.read_coord()