import numpy as np

from pywfn import _core


class Diis:
    """Pulay DIIS（Direct Inversion in the Iterative Subspace）外推器"""

    def __init__(self, max_len: int = 8) -> None:
        self.core = _core.qmcal.diis.Diis(max_len)  # type: ignore

    def reset(self) -> None:
        """清空历史向量"""
        self.core.reset()

    @staticmethod
    def commutator(fock: np.ndarray, dmat: np.ndarray, smat: np.ndarray) -> np.ndarray:
        """计算 Fock 矩阵与密度矩阵的交换子 ``[F, D]``（误差矩阵）"""
        return _core.qmcal.diis.Diis.commutator(fock, dmat, smat)

    def extrapolate(self, fock: np.ndarray, err: np.ndarray) -> np.ndarray | None:
        """用历史向量外推新的 Fock 矩阵；历史不足时返回 None"""
        return self.core.extrapolate(fock, err)
