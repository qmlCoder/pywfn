import sys
sys.path.append("D:\code\pywfn")

from pywfn.base import Mol
from pywfn.reader import LogReader


path="D:\BaiduSyncdisk\Articles\HFV\gfile\CH4.log"
reader=LogReader(path)
mol=Mol(reader=reader)
eng=reader.read_energy()
print(mol.energy)