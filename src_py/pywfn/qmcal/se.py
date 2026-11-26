from pywfn import _core
from pywfn.base import Mole


class RSE:
    def __init__(self, mole: Mole) -> None:
        self.core = _core.qmcal.se.RSE(mole.core)  # type: ignore

    def run(self) -> tuple[float, _core.base.Coefs]:  # type: ignore
        """运行半经验计算

        Returns:
            tuple[float, _core.base.Coefs]: (总能量, 分子轨道系数)，
            返回的 Coefs 为未包装的底层对象，可用 pywfn.base.Coefs(core, mole) 包装
        """
        return self.core.run()
