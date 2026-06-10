from pywfn import _core


def calc_dmat(cmat, obt_occs):
    return _core.maths.calc_dmat(cmat, obt_occs)  # type: ignore


def matrix_nhalf(mat):
    return _core.maths.matrix_nhalf(mat)  # type: ignore


def matrix_phalf(mat):
    return _core.maths.matrix_phalf(mat)  # type: ignore


def fac(n: int) -> int:
    return _core.maths.fac(n)  # type: ignore


def fac2(n: int) -> int:
    return _core.maths.fac2(n)  # type: ignore
