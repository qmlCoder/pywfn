import sys;sys.path.append('d:/code/pywfn')

from pywfn.base import Mole
from pywfn.reader import LogReader
from pywfn.fragprop import energy
from pywfn.atomprop import activity

path0=rf"d:\code\wmview_api\gfile\NHC-C4H3N3.log"
pathN=rf"d:\code\wmview_api\gfile\NHC-C4H3N3N.log"
pathP=rf"d:\code\wmview_api\gfile\NHC-C4H3N3P.log"

mol0=Mole(LogReader(path0))
molN=Mole(LogReader(pathN))
molP=Mole(LogReader(pathP))

caler=activity.Calculator(mol0)

res=caler.dirFukui(8,None,molN,molP,'mulliken')
print(res)