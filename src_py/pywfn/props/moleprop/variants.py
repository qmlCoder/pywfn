from pywfn import _core
from pywfn.base.mole import Mole


class Calculator:
    """分子变体构造器

    支持自旋分离（α/β）与 π 电子投影，生成新的分子对象。
    """

    def __init__(self, mole: Mole, atms: list[int] | None = None) -> None:
        self.mole = mole
        self.caler = _core.props.moleprop.variants.Calculator(mole.core, atms)  # type: ignore

    def mole_spin(self) -> tuple[Mole, Mole]:
        """将分子分离为 α 和 β 自旋两个独立的分子对象"""
        mole_a, mole_b = self.caler.mole_spin()
        return Mole(mole_a), Mole(mole_b)

    def mole_pi_pocv(self) -> Mole:
        """生成仅包含 π 电子成分（POCV 方法）的分子对象"""
        return Mole(self.caler.mole_pi_pocv())

    def mole_pi_mocv(self) -> Mole:
        """生成仅包含 π 电子成分（MOCV 方法）的分子对象"""
        return Mole(self.caler.mole_pi_mocv())
