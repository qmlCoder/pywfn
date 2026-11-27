"""
cub 立方体文件读取器

除了分子几何信息之外，还能读到格点的起点、步长、格点数以及格点上的数值。
"""

from pywfn import _core
from pywfn.base.atoms import Atoms
from pywfn.base.basis import Basis
from pywfn.base.coefs import Coefs


class CubReader:  # type: ignore
    """cub 立方体文件读取器"""

    def __init__(self, _core: _core.reader.CubReader):  # type: ignore
        """包装 `_core` 的读取器

        Args:
            _core: `_core.reader.CubReader` 对象
        """
        self.core = _core

    @staticmethod
    def from_path(path: str) -> "CubReader":
        """由文件路径创建读取器

        Args:
            path: .cub 文件的路径

        Returns:
            CubReader: 读取器实例
        """
        return CubReader(_core.reader.CubReader(path))  # type: ignore

    def get_atoms(self) -> "Atoms":
        """获取分子几何信息

        Returns:
            Atoms: 原子集合（`_core` 层对象）
        """
        return self.core.get_atoms()

    def get_basis(self) -> "Basis | None":
        """获取基组信息 (cub文件可能不包含基组信息)

        Returns:
            Basis | None: 基组（`_core` 层对象），文件中没有时为 None
        """
        return self.core.get_basis()

    def get_coefs(self) -> "Coefs | None":
        """获取轨道系数 (cub文件可能不包含轨道系数)

        Returns:
            Coefs | None: 轨道系数（`_core` 层对象），文件中没有时为 None
        """
        return self.core.get_coefs()

    def read_cube(self) -> "_core.props.gridprop.CubeGrid":
        """读取格点信息

        返回的是 `_core` 的立方体格点对象，可以直接用 `get()` 取坐标与数量；
        它的 `shape`/`start`/`steps` 属性同样可以通过 `.core` 访问
        （Python 层的 `props.gridprop.CubeGrid` 目前只封装了网格设置与取值，
        强行包装会丢失这三个属性）。

        Returns:
            _core.props.gridprop.CubeGrid: 立方体格点对象
        """
        return self.core.read_cube()

    def read_atms(self) -> tuple[list[int], list[float], list[list[float]]]:
        """读取原子信息

        Returns:
            tuple[list[int], list[float], list[list[float]]]:
            (原子序数列表, 原子电荷列表, 原子坐标列表)
        """
        atms, crgs, xyzs = self.core.read_atms()
        return atms, crgs, xyzs

    def read_vals(self) -> tuple[list[int], list[list[float]]]:
        """读取格点数值

        Returns:
            tuple[list[int], list[list[float]]]:
            (轨道索引列表, 每个轨道在全部格点上的数值列表)
        """
        obts, vals = self.core.read_vals()
        return obts, vals
