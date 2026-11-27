from pywfn import _core
from pywfn.base import Mole


class Calculator:
    def __init__(self, mole: Mole) -> None:
        self.mole = mole
        self.core = _core.props.orbtprop.obtval.Calculator(mole.core)  # type: ignore

    def pi_eles_pocv(self) -> list[float]:
        """使用pocv方法计算pi电子数"""
        return self.core.pi_eles_pocv()

    def pi_eles_mocv(self) -> list[float]:
        """使用mocv方法计算pi电子数"""
        return self.core.pi_eles_mocv()

    def eles(self) -> list[float]:
        """计算各原子轨道上的电子数"""
        return self.core.eles()
