from pywfn.gridprop import CubeGrid, wfnfunc
from pywfn.writer import CubWriter
from pywfn.base import Mole
from pathlib import Path
from pywfn.datas.consts import BOHR
import numpy as np
root = Path(__file__).parent


def test_cub_writer():
    mole = Mole.from_file(rf"{root}/moles/C2H4.out")
    p0, p1 = mole.border()
    cube = CubeGrid()
    cube.set_v1(p0, p1, 0.2, 4.0)
    shape, grids = cube.get()
    caler = wfnfunc.Calculator(mole)
    wfns = caler.obt_wfn(7, grids, 0)[0].reshape(1, -1)
    writer = CubWriter(rf"{root}/moles/C2H4.cub")
    writer.read_mole(mole)
    x0, y0, z0 = p0
    writer.set(
        title="pywfn生成cub文件",
        obts=[10],
        pos0=[x0-4.0, y0-4.0, z0-4.0],
        size=shape,
        step=[0.2, 0.2, 0.2],
        vals=wfns
    )
    writer.save()


def test_fch_writer():
    mole = Mole.from_file(rf"{root}/moles/C2H4.fch")


if __name__ == "__main__":
    test_cub_writer()
