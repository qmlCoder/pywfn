import sys;sys.path.append("d:/code/pywfn")

from pywfn.base import Mole
from pywfn.reader import LogReader
from pywfn.fragprop import energy


path=rf"d:\gfile\静电相互作用\fh_fh.log"

mol=Mole(LogReader(path))

caler=energy.Calculator(mol)
res=caler.EIEBA([1,2],[3,4])
print(res)

print(sum(res))