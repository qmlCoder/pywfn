
import numpy as np
from pywfn.base import Mole
from pywfn.atomprop import charge, activity, direction
from pathlib import Path


root = Path(__file__).parent


def test_charge():
    path = f"{root}/moles/C2H4.fch"
    mole = Mole.from_file(path)
    caler = charge.Calculator(mole)
    res = caler.mulliken()
    print(res)
    res = caler.lowdin()
    print(res)
    res = caler.hirshfeld()
    ref = np.array([-1.14231901, -1.19689472, -0.10210515, -0.23701916, -1.04734702, -
                   0.75745533, 0.66749594, 0.69583911, 0.88428708, 0.85552831, 0.58581315, 0.79417679])
    print(res)
    res = caler.pi_pocv()
    print(res)
    res = caler.pi_mocv()
    print(res)

def test_activity():
    # 测试福井函数
    mole_0 = Mole.from_file(f"{root}/moles/C6H6.fch")
    mole_n = Mole.from_file(f"{root}/moles/C6H6_N.fch")
    mole_p = Mole.from_file(f"{root}/moles/C6H6_P.fch")
    caler=activity.Calculator(mole_0)
    res=caler.fukui(mole_n,mole_p)
    print("fukui:\n",res)
    res=caler.fukui_pi(mole_n,mole_p)
    print("fukui_pi:\n",res)
    res=caler.fukui_dir(0,[0,0,1],mole_n,mole_p)
    print("fukui_dir:\n",res)
    # 测试自由价
    res=caler.freev()
    print("freev:\n",res)
    res=caler.freev_pi()
    print("freev_pi:\n",res)
    res=caler.freev_dir(0,[0,0,1],[1,6,5])
    print("freev_dir:\n",res) 

def test_direction():
    mole=Mole.from_file(rf"{root}/moles/test_dir.fch")
    caler=direction.Calculator(mole)
    res=caler.normal_vector(0)
    print("法向量:\n",res)
    res=caler.reactions(7)
    print("反应方向:\n",res)
    res=caler.LCS(5)
    print("局部坐标系:\n",res)

if __name__=="__main__":
    # test_activity()
    test_direction()
