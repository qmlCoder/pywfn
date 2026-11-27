"""
定义gjf的读取器，感觉可能也不会用得到
"""

from pywfn import _core
from pywfn.base.atoms import Atoms
from pywfn.base.basis import Basis
from pywfn.base.coefs import Coefs
from pywfn.reader import GauJob, _wrap_gau_job


class GjfReader:  # type: ignore
    """Gaussian 输入文件 (.gjf) 读取器"""

    def __init__(self, _core: _core.reader.GjfReader):  # type: ignore
        """包装 `_core` 的读取器

        Args:
            _core: `_core.reader.GjfReader` 对象
        """
        self.core = _core

    @staticmethod
    def from_path(path: str) -> "GjfReader":
        """由文件路径创建读取器

        Args:
            path: .gjf 文件的路径

        Returns:
            GjfReader: 读取器实例
        """
        return GjfReader(_core.reader.GjfReader(path))  # type: ignore

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

    def read_geome(self) -> tuple[list[int], list[list[float]]]:
        """读取几何结构

        Returns:
            tuple[list[int], list[list[float]]]: (原子序数列表, 原子坐标列表)
        """
        atms, xyzs = self.core.read_geome()
        return atms, xyzs

    def read_charge(self) -> tuple[int, int]:
        """读取电荷和自旋多重度

        Returns:
            tuple[int, int]: (电荷, 自旋多重度)
        """
        charge, multip = self.core.read_charge()
        return charge, multip

    def read_job(self) -> str:
        """读取作业关键词

        Returns:
            str: 以 # 开头的路由行内容
        """
        return self.core.read_job()

    def read_gau_job(self) -> GauJob | None:
        """解析作业关键词，得到任务类型、泛函与基组

        Returns:
            GauJob | None: 解析成功时返回任务信息，否则返回 None
        """
        core_job = self.core.read_gau_job()
        return None if core_job is None else _wrap_gau_job(core_job)
