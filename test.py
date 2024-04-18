from pywfn.base import Mol
from pywfn.reader import LogReader
from pywfn.bondprop import mayer

# pathn="D:\BaiduSyncdisk\Articles\HFV\gfile\CH4-.log"
# path0="D:\BaiduSyncdisk\Articles\HFV\gfile\CH4.log"
# pathp="D:\BaiduSyncdisk\Articles\HFV\gfile\CH4+.log"

root="D:\BaiduSyncdisk\Articles\HFV\gfile\M4\M4_wfn"
paths=[f'{root}{tail}.log' for tail in ['-','','+']]

mols=[Mol(reader=LogReader(path)) for path in paths]

cals=[mayer.Calculator(mol) for mol in mols]

bonds={
    3:[4,2,25],
    27:[26,28,45]
}

for c,cal in enumerate(cals):
    print(['-1价','0价','+1价'][c])
    for a1 in bonds.keys():
        orders=[]
        for a2 in bonds[a1]:
            cal.bond=a1,a2
            order=cal.calculate()
            orders.append(order)
            print(f'{a1}->{a2}:{order:.4f}')
        print(f'键级之和:{sum(orders):.4f}')