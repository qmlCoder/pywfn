"""
此脚本用来读取fchk文件
fchk文件中有哪些属性是可以用到的？
"""

import numpy as np

from pywfn.base.basis import Basis
from pywfn.base.coefs import Coefs
from pywfn.base.atoms import Atoms, Stm
from pywfn import _core
from pywfn.reader import _wrap_stm


class FchReader:  # type: ignore
    """Gaussian 格式化检查点 (.fch/.fchk) 文件读取器"""

    def __init__(self, _core: _core.reader.FchReader): # type: ignore
        """包装 `_core` 的读取器

        Args:
            _core: `_core.reader.FchReader` 对象
        """
        self.core = _core
    
    @staticmethod
    def from_path(path: str) -> "FchReader":
        """由文件路径创建读取器

        Args:
            path: .fch/.fchk 文件的路径

        Returns:
            FchReader: 读取器实例
        """
        return FchReader(_core.reader.FchReader(path)) # type: ignore

    def get_atoms(self) -> "Atoms":
        """获取原子集合

        Returns:
            Atoms: 原子集合（`_core` 层对象）
        """
        return self.core.get_atoms()

    def get_basis(self) -> "Basis | None":
        """获取基组

        Returns:
            Basis | None: 基组（`_core` 层对象），文件中没有时为 None
        """
        return self.core.get_basis()

    def get_coefs(self) -> "Coefs | None":
        """获取轨道系数

        Returns:
            Coefs | None: 轨道系数（`_core` 层对象），文件中没有时为 None
        """
        return self.core.get_coefs()

    def read_atms(self) -> list[int]:
        """读取原子序数

        Returns:
            list[int]: 每个原子的原子序数（带符号，用于标记赝势等特殊原子）
        """
        return self.core.read_atms()

    def read_crgs(self) -> list[float]:
        """读取原子核电荷

        Returns:
            list[float]: 每个原子的核电荷（无 ECP 时即原子序数）
        """
        return self.core.read_crgs()

    def read_xyzs(self) -> list[list[float]]:
        """读取笛卡尔坐标

        Returns:
            list[list[float]]: 每个原子的 [x, y, z] 坐标
        """
        return self.core.read_xyzs()

    def read_stms(self) -> list[Stm] | None:
        """读取局部坐标系

        Returns:
            list[Stm] | None: 每个原子的局域坐标系，文件中没有时为 None
        """
        core_stms = self.core.read_stms()
        return None if core_stms is None else [_wrap_stm(stm) for stm in core_stms]

    def read_coefs(
        self,
    ) -> tuple[
        list[int], list[int], list[str], list[float], list[float], np.ndarray
    ]:
        """读取轨道系数及相关的轨道信息

        Returns:
            tuple[list[int], list[int], list[str], list[float], list[float], np.ndarray]:
            (基函数所属原子索引, 基函数所属壳层索引, 基函数符号,
            轨道能量列表, 轨道占据数列表, 轨道系数矩阵)
        """
        return self.core.read_coefs()
