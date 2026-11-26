import numpy as np
from pywfn import _core
from pywfn.base import Mole


class RHF:
    def __init__(self, mole: Mole) -> None:
        self.core = _core.qmcal.hf.RHF(mole.core)  # type: ignore

    def run(self, emat: np.ndarray | None = None) -> None:
        """运行Hartree-Fock计算

        Args:
            emat: 可选，双电子积分矩阵（四维数组）。不提供时 Rust 侧自行计算
        """
        self.core.run(emat)
