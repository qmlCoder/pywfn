from pywfn import _core
from pywfn.base.mole import Mole


class Calculator:
    """拓扑分析计算器

    从若干初猜点出发寻找分子电子密度（或 π 密度）的临界点，用于 AIM 等拓扑分析。
    """

    def __init__(self, mole: Mole) -> None:
        self.mole = mole
        self.caler = _core.props.gridprop.topolog.Calculator(mole.core)  # type: ignore

    def get_grids(self) -> list[list[float]]:
        """获取分子中的临界点初猜（原子中心 + 键中心）"""
        return self.caler.get_grids()

    def analysis(self, func: str = "mol_rho_dm") -> list[tuple[list[float], float, int]]:
        """执行拓扑分析，返回全部临界点

        Args:
            func: 密度函数名，可选 "mol_rho_dm" / "pi_rho_pocv" / "pi_rho_mocv"

        Returns:
            每个临界点的 ``(位置 [x,y,z], 密度值, 临界点类型)``；
            类型为 Hessian 正、负特征值个数之差（3=核临界点，-3=笼临界点等）。
        """
        return self.caler.analysis(func)

    def get_point(
        self, grid: list[float], func: str = "mol_rho_dm"
    ) -> tuple[list[float], float, int]:
        """从指定格点出发寻找临界点

        Args:
            grid: 初始格点坐标
            func: 密度函数名，可选 "mol_rho_dm" / "pi_rho_pocv" / "pi_rho_mocv"

        Returns:
            (位置, 密度值, 临界点类型)
        """
        return self.caler.get_point(grid, func)

    def get_func(self, func: str = "mol_rho_dm"):
        """获取指定名称的密度函数对象

        Args:
            func: 密度函数名，可选 "mol_rho_dm" / "pi_rho_pocv" / "pi_rho_mocv"

        Returns:
            可调用对象，输入格点返回函数值
        """
        return self.caler.get_func(func)
