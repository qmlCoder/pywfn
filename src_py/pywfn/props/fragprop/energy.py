from pywfn import _core
from pywfn.base.mole import Mole


class Calculator:
    def __init__(self, mole: Mole) -> None:
        self.mole = mole
        self.core = _core.props.fragprop.energy.Calculator(mole.core)  # type: ignore

    def EIEBA(self, fragA: list[int], fragB: list[int]) -> tuple[float, float, float]:
        """计算两个碎片之间的相互作用能

        Args:
            fragA: 碎片 A 的原子索引列表
            fragB: 碎片 B 的原子索引列表

        Returns:
            tuple[float, float, float]: (静电相互作用能, 交换排斥能, 轨道相互作用能)
        """
        return self.core.EIEBA(fragA, fragB)
