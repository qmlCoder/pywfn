"""与 rswfn 对齐后新增/补齐接口的测试。

覆盖：utils、maths、datas、space、editor、base、props、reader、rsxcs。
只做轻量计算，不跑 SCF。
"""

from pathlib import Path

import numpy as np
import pytest

import pywfn
from pywfn import _core
from pywfn.base import Atoms, Mole, Vec3

root = Path(__file__).parent.parent


def test_version():
    assert isinstance(pywfn.rswfn_version(), str)
    assert pywfn.rswfn_version() == _core.utils.rswfn_version()


def test_utils():
    from pywfn import utils

    assert utils.file_name(r"E:\a\b.txt") == "b"
    assert utils.file_type(r"E:\a\b.txt") == "txt"
    assert utils.fort_form(1.23456, 3) == "1.235E+00"
    assert utils.all_same([1.0, 2.0], [1.0, 2.0]) is True
    assert utils.all_same([1.0, 2.0], [1.0, 3.0]) is False

    mat = np.array([[1.0, 2.0], [2.0, 4.0]])
    assert utils.symm_matrix(mat) is True
    assert utils.posi_matrix(np.eye(2)) is True
    assert utils.posi_matrix(np.array([[1.0, 2.0], [2.0, -1.0]])) is False
    assert utils.has_nan(mat) is False
    assert utils.same_matrix(mat, mat.copy()) is True
    utils.show_matrix(mat)
    utils.comp_matrix(mat, mat.copy(), 1)


def test_maths_basic():
    from pywfn import maths

    assert maths.arg_max([1.0, 5.0, 3.0]) == (1, 5.0)
    assert maths.arg_min([1.0, 5.0, 3.0]) == (0, 1.0)
    assert np.allclose(maths.vec_cross([1, 0, 0], [0, 1, 0]), [0, 0, 1])
    assert np.isclose(maths.vec_len([3, 4, 0]), 5.0)
    assert np.isclose(maths.vec_dist([0, 0, 0], [3, 4, 0]), 5.0)
    assert np.allclose(maths.vec_norm([3, 4, 0]), [0.6, 0.8, 0])
    assert isinstance(maths.vec_psd([1, 2, 2]), float)
    assert np.isclose(maths.calc_dist([0, 0, 0], [3, 4, 0]), 5.0)
    assert np.allclose(maths.vec_project([1, 1, 0], [1, 0, 0]), [1, 0, 0])
    assert maths.sign(2.0, -1.0) == -2.0

    mat = np.array([[2.0, 0.0], [0.0, 4.0]])
    assert np.isclose(maths.det(mat), 8.0)
    assert np.allclose(maths.inv(mat), [[0.5, 0], [0, 0.25]])
    vals, vecs = maths.eigh(mat)
    assert np.allclose(sorted(vals), [2.0, 4.0])

    plan = maths.get_plan_by_3points(Vec3.new(0, 0, 0), Vec3.new(1, 0, 0), Vec3.new(0, 1, 0))
    assert np.allclose(plan[:3], [0, 0, 1])
    assert np.isclose(plan[3], 0.0)

    prob = maths.bolzman_dist([0.0, 1.0], 300.0)
    assert np.isclose(sum(prob), 1.0)


def test_datas():
    from pywfn import datas
    from pywfn.datas import btrans, consts, dhmo2, integ, rdens, roots, rys_root, temps

    assert np.isclose(consts.BOHR, 1.8897261254578281)
    assert integ.RHM.shape == integ.WHM.shape
    assert consts.BOHR > 0
    assert len(datas.FACS) == 6
    assert isinstance(temps.GJF, str) and len(temps.GJF) > 0

    assert len(roots.hermite(3)) == 3
    assert len(roots.legendre(3)) == 3
    assert len(roots.rys_roots(3, 0.5)) == 3
    assert len(roots.RYS_ROOTS_1) == 1
    assert len(roots.RYS_ALPS) == 50

    r1, w1 = rys_root.rys_root1(0.5)
    assert len(r1) == len(w1) == 1
    assert len(rys_root.rys_roots(4, 0.5)) == 4
    assert len(rys_root.POLY_LARGEX_RT) == 496

    assert len(dhmo2.get_paras(6, 6, 2, 2)) == 59
    assert rdens.get_rad_dens_v1(6, 0.5) >= 0

    assert btrans.lmn2sym([1, 0, 0]) == "PX"
    assert btrans.sym2lmn("PX") == [1, 0, 0]
    assert btrans.sym2ang("PX") == 1
    assert btrans.sym2stp("S") == 0
    assert len(btrans.ang2lmns(2)) == 6
    assert len(btrans.ang2syms(2, "car")) == 6
    assert btrans.get_sph_syms() == btrans.get_sph_syms_no_cache()
    assert len(btrans.len2syms(1)) == 1

    from pywfn.datas import elems

    assert elems.atm2sym_no_cache(6) == "C"
    assert elems.get_by_idx_no_cache(6).sym == "C"
    assert np.isclose(elems.zeff_no_cache(6, 1), elems.zeff(6, 1))


def test_space_and_editor():
    from pywfn import editor, space

    grids = space.get_grids(2, 2, 2)
    assert grids.shape == (8, 3)

    pts = np.array([[1.0, 0.0, 0.0]])
    res = editor.rotate_points(pts, [0, 0, 0], [0, 0, 1], np.pi / 2)
    assert np.isclose(np.linalg.norm(res[0]), 1.0)
    assert np.isclose(res[0][2], 0.0)
    assert np.allclose(np.abs(res[0][:2]), [0.0, 1.0], atol=1e-8)


def test_base_extras():
    from pywfn.base import Ato, Bas

    ato = Ato(_core.base.Ato(0, 0, 1, "PX", 0, "sph", [1.0], [1.0]))
    assert ato.ang == 1
    assert ato.lmn() == [1, 0, 0]
    bas = Bas(_core.base.Bas(0, 0, 1, 1.0, 1.0))
    assert bas.alp == 1.0

    atoms = _core.base.Atoms([6, 1], [[0, 0, 0], [1, 0, 0]], [0.0, 0.0])
    assert atoms.len() == 2
    assert np.isclose(atoms.dmat()[0, 1], 1.0)
    assert np.isclose(atoms.cent([0, 1]).x, 0.5)
    assert len(atoms.idxs()) == 2
    assert atoms.nuc_eng() > 0.0

    mole = Mole.from_file(f"{root}/moles/C6H6.fch")
    assert mole.formula() != ""
    assert mole.ncrg() > 0
    assert len(mole.rings()) >= 1
    assert mole.equal(mole.clone()) is True


def test_props_extras():
    from pywfn.props.atomprop import charge
    from pywfn.props.bondprop import order
    from pywfn.props.gridprop import GridVals, dftgrid
    from pywfn.props.orbtprop import mocv

    gv = GridVals(1.0, [0.0, 1.0, 0.0], [[0.0] * 3] * 3)
    assert gv.val0 == 1.0
    assert np.allclose(gv.val1, [0.0, 1.0, 0.0])
    assert len(gv.val2s()) == 9

    mole = Mole.from_file(f"{root}/moles/C6H6.fch")
    caler = dftgrid.Calculator(mole)
    grids, weits, meta = caler.atm_grids_meta(0, 10, 14, True)
    assert grids.shape[0] == len(weits) == len(meta)

    dmat = mole.get_dmat("sph")
    smat = mole.get_smat("sph")
    omat = order.order(dmat, smat, "mayer")
    assert omat.shape == dmat.shape
    uplo = mole.basis.shls.atm_uplo("sph")
    eles = charge.qmat2eles(dmat, uplo)
    assert len(eles) == mole.atoms.len()

    assert mocv.OrderDType.PI != mocv.OrderDType.SI
    dest = mocv.MocvDest.OrderDecom(mocv.OrderDType.PI)
    assert isinstance(mocv.mocv_atos(mole, dest), list)


def test_resp():
    from pywfn.props.atomprop import resp

    mole = Mole.from_file(f"{root}/moles/C2H4.out")
    caler = resp.Calculator(mole)
    res = caler.resp()
    assert len(res.charges) == mole.atoms.len()
    assert res.n_points > 0
    assert np.all(np.isfinite(res.charges))
    assert abs(res.charges.sum()) < 0.5

    opts = resp.RespOptions.default()
    assert opts.pt_density > 0
    res2 = caler.resp_with(opts)
    assert np.allclose(res2.charges, res.charges, atol=1e-5)

    pts = caler.mk_points(opts)
    assert pts.shape[1] == 3
    vals = caler.esp_values(pts[:5])
    assert vals.shape == (5,)


def test_reader_extras():
    from pywfn.reader import FileType, charge_number, get_fileType, read_lines

    assert get_fileType(r"E:\a\b.fchk") == FileType.Fch
    assert get_fileType(r"E:\a\b.out") == FileType.Log
    assert charge_number(0, 1, 16) == (8, 8)
    lines = read_lines(f"{root}/moles/C6H6.fch")
    assert len(lines) > 5

    mole = Mole.from_file(f"{root}/moles/C6H6.fch")
    assert mole.charge == 0


def test_rsxcs():
    from pywfn import rsxcs

    names = rsxcs.supported_names()
    assert "lda_x" in names
    func = _core.rsxcs.XcFunc.from_name("b3lyp")
    assert func is not None
    assert len(rsxcs.supported()) > 0


def test_march():
    from pywfn.maths import march

    n = 5
    shape = [n, n, n]
    xs = np.linspace(-1, 1, n)
    grids = np.array([[x, y, z] for x in xs for y in xs for z in xs])
    values = (grids**2).sum(axis=1)
    verts, types = march.march_cube_v1(shape, grids, values, 0.5)
    assert len(verts) == len(types)
    verts2, types2 = march.surface_nets_v1(shape, grids, values, 0.5)
    assert len(verts2) == len(types2)


if __name__ == "__main__":
    test_utils()
    test_maths_basic()
    test_datas()
    test_space_and_editor()
    test_base_extras()
    test_props_extras()
    test_reader_extras()
    test_rsxcs()
    test_march()
