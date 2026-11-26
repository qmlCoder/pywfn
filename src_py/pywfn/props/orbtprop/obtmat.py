import numpy as np

from pywfn import _core
from pywfn.base import Mole, Stm


class Calculator:
    def __init__(self, mole: Mole, atms: list[int] | None = None) -> None:
        self.mole = mole
        self.core = _core.props.orbtprop.obtmat.Calculator(mole.core, atms)  # type: ignore

    def pocv(
        self, dirs: dict[int, list[float]], keep_other_atm: bool, keep_other_sym: bool
    ) -> np.ndarray:
        """计算投影分子轨道的系数矩阵

        Args:
            dirs: 方向字典, {原子索引: 方向向量}
            keep_other_atm: 是否保留其他原子
            keep_other_sym: 是否保留其他对称性

        Returns:
            np.ndarray: 投影后的系数矩阵
        """
        return self.core.pocv(dirs, keep_other_atm, keep_other_sym)

    def pi_pocv(self) -> tuple[dict[int, list[float]], np.ndarray]:
        """计算pi轨道的投影系数矩阵

        Returns:
            tuple: (每个原子的法向量方向, 投影后的系数矩阵)
        """
        return self.core.pi_pocv()

    def mocv(self, stms: dict[int, "Stm"], atos: list[int]) -> np.ndarray:
        """使用MOCV方法分解轨道系数矩阵

        Args:
            stms: 原子局部坐标系字典
            atos: 原子轨道索引

        Returns:
            np.ndarray: 分解后的系数矩阵
        """
        stms = {atm: stm.core for atm, stm in stms.items()}
        return self.core.mocv(stms, atos)

    def pi_mocv(
        self, dest: str = "PiProperty"
    ) -> tuple[dict[int, Stm], list[int], np.ndarray]:
        """计算pi轨道的MOCV分解

        Args:
            dest: 分解目标, 可选 "PiProperty"、"OrderDecom-SI"、"OrderDecom-PI"、
                "OrderDecom-DE"、"OrderDecom-PH"、"RoundOrder"

        Returns:
            tuple: (局部坐标系, 保留的原子轨道索引, 分解后的系数矩阵)
        """
        stms, atos, cmat = self.core.pi_mocv(dest)
        stms = {atm: Stm(stm) for atm, stm in stms.items()}
        return stms, atos, cmat

    def ldao(self, stms: dict[int, Stm]) -> np.ndarray:
        """使用LDAO方法分解轨道系数矩阵

        Args:
            stms: 原子局部坐标系字典

        Returns:
            np.ndarray: 分解后的系数矩阵
        """
        stms = {atm: stm.core for atm, stm in stms.items()}
        return self.core.ldao(stms)


def mocv(L: np.ndarray, rcoefs: np.ndarray, keeps: np.ndarray) -> np.ndarray:
    """计算MOCV分解

    Args:
        L: 变换矩阵
        rcoefs: 系数矩阵
        keeps: 保留掩码矩阵, 与rcoefs形状相同, 值为0/1, 与rcoefs逐元素相乘

    Returns:
        np.ndarray: 分解结果
    """
    keeps = np.asarray(keeps, dtype=np.float64)
    return _core.props.orbtprop.obtmat.mocv(L, rcoefs, keeps)  # type: ignore


def mocv_atos(mole: Mole, dest: str) -> list[int]:
    """将MOCV分解结果转换为原子轨道索引

    Args:
        mole: 分子对象
        dest: PiProperty|OrderDecom-SI|OrderDecom-PI|OrderDecom-DE|RoundOrder

    Returns:
        list[int]: 原子轨道索引
    """
    return _core.props.orbtprop.obtmat.mocv_atos(mole.core, dest)  # type: ignore
