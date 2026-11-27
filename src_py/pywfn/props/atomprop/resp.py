"""
RESP 电荷

对应 rswfn 的 `props::atomprop::resp` 模块，实现约束静电势拟合。
"""

from enum import IntEnum

import numpy as np

from pywfn import _core
from pywfn.base import Mole


class EspType(IntEnum):
    """静电势的类型"""

    Total = 0  # 总静电势（核 + 电子）
    Electronic = 1  # 仅电子贡献


_CORE_ESP = {
    EspType.Total: _core.props.atomprop.resp.EspType.Total,  # type: ignore
    EspType.Electronic: _core.props.atomprop.resp.EspType.Electronic,  # type: ignore
}


class RespOptions:
    """RESP 拟合的参数设置

    不传参数时使用底层 rswfn 的默认参数；只传部分参数时，未指定的项沿用默认值。
    """

    def __init__(
        self,
        layer_scales: list[float] | None = None,
        pt_density: float | None = None,
        hyper_a1: float | None = None,
        hyper_a2: float | None = None,
        hyper_b: float | None = None,
        max_iter: int | None = None,
        conv: float | None = None,
        two_stage: bool | None = None,
        esp_type: EspType | None = None,
    ) -> None:
        if (
            layer_scales is None
            and pt_density is None
            and hyper_a1 is None
            and hyper_a2 is None
            and hyper_b is None
            and max_iter is None
            and conv is None
            and two_stage is None
            and esp_type is None
        ):
            self.core = _core.props.atomprop.resp.RespOptions.default()  # type: ignore
            return
        base = _core.props.atomprop.resp.RespOptions.default()  # type: ignore
        self.core = _core.props.atomprop.resp.RespOptions(  # type: ignore
            base.layer_scales if layer_scales is None else layer_scales,
            base.pt_density if pt_density is None else pt_density,
            base.hyper_a1 if hyper_a1 is None else hyper_a1,
            base.hyper_a2 if hyper_a2 is None else hyper_a2,
            base.hyper_b if hyper_b is None else hyper_b,
            base.max_iter if max_iter is None else max_iter,
            base.conv if conv is None else conv,
            base.two_stage if two_stage is None else two_stage,
            _CORE_ESP[EspType(esp_type)] if esp_type is not None else base.esp_type,
        )

    @staticmethod
    def default() -> "RespOptions":
        """使用 rswfn 的默认参数"""
        opts = RespOptions.__new__(RespOptions)
        opts.core = _core.props.atomprop.resp.RespOptions.default()  # type: ignore
        return opts


    @property
    def layer_scales(self) -> list[float]:
        """各层格点密度的缩放系数"""
        return self.core.layer_scales

    @layer_scales.setter
    def layer_scales(self, val: list[float]) -> None:
        self.core.layer_scales = val

    @property
    def pt_density(self) -> float:
        """格点密度"""
        return self.core.pt_density

    @pt_density.setter
    def pt_density(self, val: float) -> None:
        self.core.pt_density = val

    @property
    def hyper_a1(self) -> float:
        """第一层双曲惩罚系数"""
        return self.core.hyper_a1

    @hyper_a1.setter
    def hyper_a1(self, val: float) -> None:
        self.core.hyper_a1 = val

    @property
    def hyper_a2(self) -> float:
        """第二层双曲惩罚系数"""
        return self.core.hyper_a2

    @hyper_a2.setter
    def hyper_a2(self, val: float) -> None:
        self.core.hyper_a2 = val

    @property
    def hyper_b(self) -> float:
        """双曲惩罚的偏移量"""
        return self.core.hyper_b

    @hyper_b.setter
    def hyper_b(self, val: float) -> None:
        self.core.hyper_b = val

    @property
    def max_iter(self) -> int:
        """最大迭代次数"""
        return self.core.max_iter

    @max_iter.setter
    def max_iter(self, val: int) -> None:
        self.core.max_iter = val

    @property
    def conv(self) -> float:
        """收敛阈值"""
        return self.core.conv

    @conv.setter
    def conv(self, val: float) -> None:
        self.core.conv = val

    @property
    def two_stage(self) -> bool:
        """是否使用两阶段拟合"""
        return self.core.two_stage

    @two_stage.setter
    def two_stage(self, val: bool) -> None:
        self.core.two_stage = val

    @property
    def esp_type(self) -> EspType:
        """静电势类型"""
        core = self.core.esp_type
        if core == _CORE_ESP[EspType.Electronic]:
            return EspType.Electronic
        return EspType.Total

    @esp_type.setter
    def esp_type(self, val: EspType) -> None:
        self.core.esp_type = _CORE_ESP[EspType(val)]


class RespResult:
    """RESP 拟合结果"""

    def __init__(self, core) -> None:
        self.core = core

    @property
    def charges(self) -> np.ndarray:
        """拟合得到的电荷"""
        return np.array(self.core.charges)

    @property
    def stage1(self) -> np.ndarray:
        """第一阶段得到的电荷"""
        return np.array(self.core.stage1)

    @property
    def n_points(self) -> int:
        """拟合使用的格点数量"""
        return self.core.n_points

    @property
    def rmse(self) -> float:
        """均方根误差"""
        return self.core.rmse

    @property
    def rrmse(self) -> float:
        """相对均方根误差"""
        return self.core.rrmse

    @property
    def stage1_iter(self) -> int:
        """第一阶段的迭代次数"""
        return self.core.stage1_iter

    @property
    def stage2_iter(self) -> int:
        """第二阶段的迭代次数"""
        return self.core.stage2_iter

    @property
    def stage2(self) -> bool:
        """是否进行了第二阶段拟合"""
        return self.core.stage2

    @property
    def converged(self) -> bool:
        """是否收敛"""
        return self.core.converged


class Calculator:
    """RESP 电荷计算器"""

    def __init__(self, mole: Mole, atms: list[int] | None = None) -> None:
        self.mole = mole
        self.caler = _core.props.atomprop.resp.Calculator(mole.core, atms)  # type: ignore

    def resp(self) -> RespResult:
        """使用默认参数计算 RESP 电荷"""
        return RespResult(self.caler.resp())

    def resp_with(self, opts: RespOptions) -> RespResult:
        """使用指定参数计算 RESP 电荷

        Args:
            opts: 拟合参数
        """
        return RespResult(self.caler.resp_with(opts.core))

    def mk_points(self, opts: RespOptions) -> np.ndarray:
        """生成用于拟合的格点

        Args:
            opts: 拟合参数
        """
        return np.array(self.caler.mk_points(opts.core))

    def esp_values(self, pts: np.ndarray, esp_type: EspType = EspType.Total) -> np.ndarray:
        """计算指定格点上的静电势

        Args:
            pts: 格点坐标
            esp_type: 静电势类型
        """
        return np.array(self.caler.esp_values(pts, _CORE_ESP[EspType(esp_type)]))
