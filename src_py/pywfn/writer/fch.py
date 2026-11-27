"""
将分子对象保存为fch文件
"""

from pywfn import _core
from pywfn.base.atoms import Atoms
from pywfn.base.basis import Basis
from pywfn.base.coefs import Coefs
from pywfn.base.mole import Mole


class FchWriter:
    """Gaussian 格式化检查点 (.fch/.fchk) 文件写出器"""

    def __init__(self) -> None:
        """fch文件写出器"""
        self.core = _core.writer.FchWriter()  # type: ignore

    @staticmethod
    def from_mole(mole: Mole) -> "FchWriter":
        """从分子对象构建写出器"""
        writer = FchWriter.__new__(FchWriter)
        writer.core = _core.writer.FchWriter.from_mole(mole.core)  # type: ignore
        return writer

    def _mole(self) -> Mole:
        """用当前写出器里的原子、基组、系数拼一个分子对象

        Python 层的 `Atoms`/`Basis`/`Coefs` 都要绑定一个分子，
        因此这里按需重建一个（不做任何文件读写）。

        Returns:
            Mole: 与写出器内容一致的分子对象
        """
        core_mole = _core.base.Mole.build(  # type: ignore
            self.core.atoms, self.core.basis, self.core.coefs
        )
        return Mole(core_mole)

    @property
    def charge(self) -> int:
        """分子总电荷"""
        return self.core.charge

    @charge.setter
    def charge(self, charge: int) -> None:
        self.core.charge = charge

    @property
    def multip(self) -> int:
        """自旋多重度"""
        return self.core.multip

    @multip.setter
    def multip(self, multip: int) -> None:
        self.core.multip = multip

    @property
    def energy(self) -> float:
        """分子能量"""
        return self.core.energy

    @energy.setter
    def energy(self, energy: float) -> None:
        self.core.energy = energy

    @property
    def atoms(self) -> Atoms:
        """原子集合"""
        return Atoms(self.core.atoms, self._mole())

    @atoms.setter
    def atoms(self, atoms: Atoms) -> None:
        self.core.atoms = atoms.core

    @property
    def basis(self) -> Basis:
        """基组信息"""
        return Basis(self.core.basis, self._mole())

    @basis.setter
    def basis(self, basis: Basis) -> None:
        self.core.basis = basis.core

    @property
    def coefs(self) -> Coefs:
        """轨道系数信息"""
        return Coefs(self.core.coefs, self._mole())

    @coefs.setter
    def coefs(self, coefs: Coefs) -> None:
        self.core.coefs = coefs.core

    def build(self) -> str:
        """生成fch文件内容"""
        return self.core.build()

    def save(self, path: str) -> None:
        """保存fch文件"""
        self.core.save(path)

    def read_mole(self, mole: Mole) -> None:
        """从分子对象中读取数据"""
        self.core.read_mole(mole.core)
