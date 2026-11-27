"""
交换关联泛函（XC）计算

对应 rswfn 的 `rsxcs` 模块，提供密度泛函的各个组成成分（LDA/GGA、交换/相关）、
泛函定义，以及它们在 RKS/UKS 下的能量密度与势的求值。

- 成分选择：`XcComponent.LdaX`、`XcComponent.B88` 等
- 泛函定义：`XcFunc.B3LYP`、`XcFunc.from_name("b3lyp")` 等
- 求值：`eval_rks` / `eval_uks`（单个成分）、`RsXcs`（整个泛函）
"""

from pywfn import _core


class XcComponent:
    """交换关联泛函的组成成分

    每个成分对应一种具体的密度泛函近似，例如 LDA 交换（LdaX）、
    相关（Vwn5 / Vwn3 / Lyp / P86 / PbeC）、GGA 交换（B88 / PbeX）等。
    直接使用类属性选择成分，例如 `XcComponent.B88`。
    """

    LdaX = _core.rsxcs.XcComponent.LdaX  # type: ignore # LDA 交换
    Vwn5 = _core.rsxcs.XcComponent.Vwn5  # type: ignore # VWN5 相关（LDA）
    Vwn3 = _core.rsxcs.XcComponent.Vwn3  # type: ignore # VWN3 相关（LDA）
    B88 = _core.rsxcs.XcComponent.B88  # type: ignore # B88 交换（GGA）
    Lyp = _core.rsxcs.XcComponent.Lyp  # type: ignore # LYP 相关（GGA）
    P86 = _core.rsxcs.XcComponent.P86  # type: ignore # P86 相关（GGA）
    PbeX = _core.rsxcs.XcComponent.PbeX  # type: ignore # PBE 交换（GGA）
    PbeC = _core.rsxcs.XcComponent.PbeC  # type: ignore # PBE 相关（GGA）

    ALL = _core.rsxcs.XcComponent.ALL  # type: ignore # 所有内置成分

    @staticmethod
    def is_gga(comp: "_core.rsxcs.XcComponent") -> bool:
        """该成分是否为广义梯度近似（GGA）

        Args:
            comp: 泛函成分，如 `XcComponent.B88`

        Returns:
            bool: GGA 成分返回 True，纯 LDA 成分返回 False
        """
        return comp.is_gga()


class XcFunc:
    """交换关联泛函

    内置 Svwn5、BLYP、BP86、PBE、PBE0、B3LYP 六种泛函，也支持通过
    `XcFunc.External(id)` 引用外部注册的泛函。
    直接使用类属性选择泛函，例如 `XcFunc.B3LYP`。
    """

    Svwn5 = _core.rsxcs.XcFunc.Svwn5  # type: ignore # Svwn5（LDA）
    BLYP = _core.rsxcs.XcFunc.BLYP  # type: ignore # BLYP（GGA）
    BP86 = _core.rsxcs.XcFunc.BP86  # type: ignore # BP86（GGA）
    PBE = _core.rsxcs.XcFunc.PBE  # type: ignore # PBE（GGA）
    PBE0 = _core.rsxcs.XcFunc.PBE0  # type: ignore # PBE0（杂化）
    B3LYP = _core.rsxcs.XcFunc.B3LYP  # type: ignore # B3LYP（杂化）

    ALL = _core.rsxcs.XcFunc.ALL  # type: ignore # 所有内置泛函

    @staticmethod
    def from_name(name: str) -> "_core.rsxcs.XcFunc | None":
        """按名称查找泛函（不区分大小写）

        Args:
            name: 泛函名称，如 "b3lyp"

        Returns:
            _core.rsxcs.XcFunc | None: 找到时返回泛函对象，否则返回 None
        """
        return _core.rsxcs.XcFunc.from_name(name)  # type: ignore

    @staticmethod
    def External(id: int) -> "_core.rsxcs.XcFunc":
        """引用外部注册的泛函

        Args:
            id: 外部泛函编号

        Returns:
            _core.rsxcs.XcFunc: 对应的泛函对象
        """
        return _core.rsxcs.XcFunc.External(id)  # type: ignore


class RksEval:
    """RKS（闭壳层）泛函求值结果

    保存各格点上的能量密度及其对密度、梯度平方的导数。

    Args:
        core: 底层 `_core` 的 RksEval 对象
    """

    def __init__(self, core: "_core.rsxcs.RksEval") -> None:
        self.core = core

    @property
    def exc(self) -> list[float]:
        """能量密度"""
        return self.core.exc

    @property
    def vrho(self) -> list[float]:
        """能量密度对密度 rho 的导数"""
        return self.core.vrho

    @property
    def vsigma(self) -> list[float]:
        """能量密度对梯度平方 sigma 的导数"""
        return self.core.vsigma

    def __repr__(self) -> str:
        return self.core.__repr__()


class UksEval:
    """UKS（开壳层）泛函求值结果

    保存各格点上的能量密度及其对 α/β 密度、梯度平方的导数。

    Args:
        core: 底层 `_core` 的 UksEval 对象
    """

    def __init__(self, core: "_core.rsxcs.UksEval") -> None:
        self.core = core

    @property
    def exc(self) -> list[float]:
        """能量密度"""
        return self.core.exc

    @property
    def vrho_a(self) -> list[float]:
        """能量密度对 α 密度 rho_a 的导数"""
        return self.core.vrho_a

    @property
    def vrho_b(self) -> list[float]:
        """能量密度对 β 密度 rho_b 的导数"""
        return self.core.vrho_b

    @property
    def vsigma_aa(self) -> list[float]:
        """能量密度对 sigma_aa 的导数"""
        return self.core.vsigma_aa

    @property
    def vsigma_bb(self) -> list[float]:
        """能量密度对 sigma_bb 的导数"""
        return self.core.vsigma_bb

    @property
    def vsigma_ab(self) -> list[float]:
        """能量密度对 sigma_ab 的导数"""
        return self.core.vsigma_ab

    def __repr__(self) -> str:
        return self.core.__repr__()


class RsXcs:
    """交换关联泛函求值器

    按泛函定义把各成分的贡献加权求和，得到能量密度与相应的势。

    Args:
        func: 泛函对象，如 `XcFunc.B3LYP`；为 None 时使用默认泛函
    """

    def __init__(self, func: "_core.rsxcs.XcFunc | None" = None) -> None:
        if func is None:
            self.core = _core.rsxcs.RsXcs()  # type: ignore
        else:
            self.core = _core.rsxcs.RsXcs.with_func(func)  # type: ignore

    @classmethod
    def with_func(cls, func: "_core.rsxcs.XcFunc") -> "RsXcs":
        """按指定泛函创建求值器

        Args:
            func: 泛函对象，如 `XcFunc.B3LYP`

        Returns:
            RsXcs: 泛函求值器
        """
        return cls(func)

    def func(self) -> "_core.rsxcs.XcFunc":
        """当前使用的泛函

        Returns:
            _core.rsxcs.XcFunc: 泛函对象
        """
        return self.core.func()

    def exact_x(self) -> float:
        """当前泛函的精确交换比例

        Returns:
            float: 0 到 1 之间的比例系数
        """
        return self.core.exact_x()

    def run(
        self, rho: list[float], sigma: list[float]
    ) -> tuple[list[float], list[float], list[float]]:
        """在闭壳层（RKS）下求值

        Args:
            rho: 各格点的电子密度
            sigma: 各格点的密度梯度平方

        Returns:
            tuple[list[float], list[float], list[float]]: (exc, vrho, vsigma)，
            均为与输入等长的列表
        """
        return self.core.run(rho, sigma)

    def run_spin(
        self,
        rho_a: list[float],
        rho_b: list[float],
        saa: list[float],
        sbb: list[float],
        sab: list[float],
    ) -> tuple[
        list[float],
        list[float],
        list[float],
        list[float],
        list[float],
        list[float],
    ]:
        """在开壳层（UKS）下求值

        Args:
            rho_a: 各格点的 α 电子密度
            rho_b: 各格点的 β 电子密度
            saa: 各格点的 sigma_aa
            sbb: 各格点的 sigma_bb
            sab: 各格点的 sigma_ab

        Returns:
            tuple: (exc, vrho_a, vrho_b, vsigma_aa, vsigma_bb, vsigma_ab)，
            均为与输入等长的列表
        """
        return self.core.run_spin(rho_a, rho_b, saa, sbb, sab)

    def __repr__(self) -> str:
        return self.core.__repr__()


def eval_rks(
    comp: "_core.rsxcs.XcComponent", rho: list[float], sigma: list[float]
) -> RksEval:
    """在闭壳层（RKS）下计算单个泛函成分的贡献

    Args:
        comp: 泛函成分，如 `XcComponent.LdaX`
        rho: 各格点的电子密度
        sigma: 各格点的密度梯度平方

    Returns:
        RksEval: 求值结果，含 exc、vrho、vsigma
    """
    return RksEval(_core.rsxcs.eval_rks(comp, rho, sigma))  # type: ignore


def eval_uks(
    comp: "_core.rsxcs.XcComponent",
    rho_a: list[float],
    rho_b: list[float],
    saa: list[float],
    sbb: list[float],
    sab: list[float],
) -> UksEval:
    """在开壳层（UKS）下计算单个泛函成分的贡献

    Args:
        comp: 泛函成分，如 `XcComponent.B88`
        rho_a: 各格点的 α 电子密度
        rho_b: 各格点的 β 电子密度
        saa: 各格点的 sigma_aa
        sbb: 各格点的 sigma_bb
        sab: 各格点的 sigma_ab

    Returns:
        UksEval: 求值结果，含 exc、vrho_a、vrho_b、vsigma_aa、vsigma_bb、vsigma_ab
    """
    return UksEval(_core.rsxcs.eval_uks(comp, rho_a, rho_b, saa, sbb, sab))  # type: ignore


def supported() -> list["_core.rsxcs.XcFunc"]:
    """返回所有内置泛函

    Returns:
        list[_core.rsxcs.XcFunc]: 泛函对象列表
    """
    return _core.rsxcs.supported()  # type: ignore


def supported_names() -> list[str]:
    """返回所有内置泛函的名称

    Returns:
        list[str]: 泛函名称列表，可直接传给 `XcFunc.from_name`
    """
    return _core.rsxcs.supported_names()  # type: ignore
