from pywfn.base import Mole
from pywfn.bondprop import order,direction
from pathlib import Path

root = Path(__file__).parent

def test_order():
    mole=Mole.from_file(rf"{root}/moles/C2H4.fch")
    caler=order.Calculator(mole)
    res=caler.mayer()
    print("mayer\n",res)
    res=caler.lowdin()
    print("lowdin\n",res)
    res=caler.wiberg()
    print("wiberg\n",res)
    res=caler.pi_pocv()
    print("pi_pocv\n",res)
    res=caler.pi_mocv()
    print("pi_mocv\n",res)
    res=caler.bound(0,[0,0,1],[1,2,5])
    print("bound\n",res)

def test_direction():
    mole=Mole.from_file(rf"{root}/moles/test_dir.fch")
    caler=direction.Calculator(mole)
    res=caler.verts(0,5)
    print("verts\n",res)
    

if __name__=="__main__":
    test_order()
    test_direction()