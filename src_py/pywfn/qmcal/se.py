"""
半经验计算

对应 rswfn 的 `qmcal::se` 模块，提供扩展 Hückel 方法（RSE）计算。
"""

from pywfn import _core
from pywfn.base import Coefs, Mole


class RSE:
    """RSE（扩展 Hückel）计算器

    用 Wolfsberg-Helmholz 近似构造 Hückel 哈密顿矩阵，一次对角化即可得到
    分子轨道与能量，常用于生成初猜或定性分析。

    Args:
        mole: 分子对象，需带有基组信息
    """

    def __init__(self, mole: Mole) -> None:
        self.mole = mole
        self.core = _core.qmcal.se.RSE(mole.core)  # type: ignore

    def run(self) -> tuple[float, Coefs]:
        """运行半经验（扩展 Hückel）计算

        Returns:
            tuple[float, Coefs]: (总能量, 分子轨道系数)
        """
        eng, core = self.core.run()
        return eng, self._result_coefs(core)

    def calc_rse(self, lcs: bool = False) -> tuple[float, Coefs]:
        """执行扩展 Hückel（RSE）计算

        Args:
            lcs: 是否使用局域坐标系（local coordinate system）构建基函数

        Returns:
            tuple[float, Coefs]: (总能量, 分子轨道系数)
        """
        eng, core = self.core.calc_rse(lcs)
        return eng, self._result_coefs(core)

    def _result_coefs(self, core: "_core.base.Coefs") -> Coefs:
        """把计算结果包装成带系数的 Coefs

        Coefs 包装类的 get_cmat/get_dmat 统一从所绑定的分子读取系数，
        因此先用计算结果重建一个分子对象，再从中取 Coefs。

        Args:
            core: 计算得到的 `_core` 轨道系数对象

        Returns:
            Coefs: 与结果系数绑定的轨道系数对象
        """
        res_mole = Mole.build(self.mole.atoms, self.mole.basis, Coefs(core, self.mole))
        return res_mole.coefs
