from pathlib import Path

import numpy as np
from pywfn.props.atomprop import activity, charge, direction
from pywfn.base import Mole

root = Path(__file__).parent.parent


def check_array(res, ref):
    assert res.shape == ref.shape
    assert np.allclose(res, ref, atol=1e-5)


def test_charge():
    path = f"{root}/moles/C2H4.out"
    mole = Mole.from_file(path)
    caler = charge.Calculator(mole)
    res = caler.mulliken()
    print(res)
    res = caler.lowdin()
    print(res)
    res = caler.hirshfeld()
    ref = np.array(
        [
            -1.14231901,
            -1.19689472,
            -0.10210515,
            -0.23701916,
            -1.04734702,
            -0.75745533,
            0.66749594,
            0.69583911,
            0.88428708,
            0.85552831,
            0.58581315,
            0.79417679,
        ]
    )
    print(res)
    res = caler.pi_pocv()
    print(res)
    res = caler.pi_mocv()
    print(res)
    mole = Mole.from_file(f"{root}/moles/C6H6_N.fch")
    caler = charge.Calculator(mole)
    res = caler.spin("mulliken")
    ref = np.array(
        [
            [2.99945, 3.11835, -0.11890],
            [3.26149, 2.93277, 0.32872],
            [3.26149, 2.93277, 0.32872],
            [2.99945, 3.11835, -0.11890],
            [3.26149, 2.93277, 0.32872],
            [3.26149, 2.93277, 0.32872],
            [0.49440, 0.49048, 0.00392],
            [0.49159, 0.51281, -0.02123],
            [0.49159, 0.51281, -0.02123],
            [0.49440, 0.49048, 0.00392],
            [0.49159, 0.51281, -0.02123],
            [0.49159, 0.51281, -0.02123],
        ]
    )
    print("-" * 100)
    check_array(res, ref)


def test_activity():
    # 测试福井函数
    mole_0 = Mole.from_file(f"{root}/moles/C6H6.fch")
    mole_n = Mole.from_file(f"{root}/moles/C6H6_N.fch")
    mole_p = Mole.from_file(f"{root}/moles/C6H6_P.fch")
    caler = activity.Calculator(mole_0)
    res = caler.fukui(mole_n, mole_p)
    print("fukui:\n", res)
    res = caler.fukui_pi_pocv(mole_n, mole_p)
    print("fukui_pi:\n", res)
    res = caler.fukui_dir(0, [0, 0, 1], mole_n, mole_p)
    print("fukui_dir:\n", res)
    # 测试自由价
    res = caler.freev()
    print("freev:\n", res)
    res = caler.freev_pi()
    print("freev_pi:\n", res)
    res = caler.freev_dir(0, [0, 0, 1], [1, 6, 5])
    print("freev_dir:\n", res)


def test_direction():
    mole = Mole.from_file(rf"{root}/moles/test_dir.fch")
    caler = direction.Calculator(mole)
    res = caler.normal_vector(0)
    print("法向量:\n", res)
    res = caler.reactions(7)
    print("反应方向:\n", res)
    res = caler.LCS(5)
    print("局部坐标系:\n", res)


if __name__ == "__main__":
    # test_activity()
    # test_direction()
    test_charge()
