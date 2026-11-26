from pywfn.props.gridprop import CubeGrid, wfnfunc
from pywfn.writer import CubWriter
from pywfn.base import Mole
from pathlib import Path
from pywfn.datas.consts import BOHR
import numpy as np
root = Path(__file__).parent.parent
output = root / "output"


def test_cub_writer():
    output.mkdir(exist_ok=True)
    mole = Mole.from_file(rf"{root}/moles/C2H4.out")
    p0, p1 = mole.border(4.0)
    cube = CubeGrid()
    cube.set_v1(p0, p1, 0.2, 4.0)
    shape, grids = cube.get()
    caler = wfnfunc.Calculator(mole)
    wfns = caler.obt_wfn(7, grids, 0)[0].reshape(1, -1)
    writer = CubWriter()
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
    writer.save(str(output / "C2H4.cub"))


def test_fch_writer():
    mole = Mole.from_file(rf"{root}/moles/C6H6.fch")


if __name__ == "__main__":
    test_cub_writer()
