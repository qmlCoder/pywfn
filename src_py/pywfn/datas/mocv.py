from pywfn import _core


def get_tmat(num, L):
    return _core.get_tmat(num, L)  # type: ignore


def get_tmat_sph(num, L):
    return _core.get_tmat_sph(num, L)  # type: ignore


def get_tmat_car(num, L):
    return _core.get_tmat_car(num, L)  # type: ignore


def tmat_s(L):
    return _core.tmat_s(L)  # type: ignore


def tmat_p(L):
    return _core.tmat_p(L)  # type: ignore


def tmat_d(L):
    return _core.tmat_d(L)  # type: ignore


def tmat_f(L):
    return _core.tmat_f(L)  # type: ignore


def tmat_g(L):
    return _core.tmat_g(L)  # type: ignore


def tmat_h(L):
    return _core.tmat_h(L)  # type: ignore
