import sys;sys.path.append("d:/code/pywfn")

from pywfn.base import Mole
from pywfn.reader import LogReader

from pywfn.atomprop import charge

path0=rf"d:\code\wmview_api\gfile\NHC-C4H3N3.log"
pathN=rf"d:\code\wmview_api\gfile\NHC-C4H3N3N.log"
pathP=rf"d:\code\wmview_api\gfile\NHC-C4H3N3P.log"

mol0=Mole(LogReader(path0))
molN=Mole(LogReader(pathN))
molP=Mole(LogReader(pathP))

mols=[mol0,molN,molP]

for mol in mols:
    caler=charge.Calculator(mol)
    res=caler.spin(ctype='mulliken')

    print(res)
