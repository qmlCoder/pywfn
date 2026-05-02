from pywfn import _core


class Stm:
    """原子局部坐标系"""

    def __init__(self, _core: _core.base.Stm) -> None:  # type: ignore
        self.core = _core

    def __repr__(self) -> str:
        return f"pywfn: {self.core.__repr__()}"


class Atom:
    """
    原子对象
    """

    def __init__(self, _core: _core.base.Atom) -> None:  # type: ignore
        self.core = _core

    def rad(self):
        return self.core.rad()

    def sym(self):
        return self.core.sym()


class Atoms:  # type: ignore

    def __init__(self, core: _core.base.Atoms) -> None:  # type: ignore
        self.core = core

    def len(self):
        return self.core.len()

    def dist(self, i: int, j: int):
        return self.core.dist(i, j)

    def syms(self):
        return self.core.syms()

    def xyzs(self):
        return self.core.xyzs()

    def __repr__(self) -> str:
        return self.core.__repr__()
