//! 交换关联泛函（XC）计算模块。
//!
//! 对应 `rswfn::rsxcs`，提供密度泛函的各成分（LDA/GGA、交换/相关）
//! 及其在 RKS/UKS 下的能量密度与势。

#![allow(non_snake_case)]

use pyo3::prelude::*;
use rswfn;

/// 交换关联泛函的组成成分
///
/// 每个成分对应一种具体的密度泛函近似，例如 LDA 交换（`LdaX`）、
/// 相关（`Vwn5` / `Lyp` / `PbeC`）、GGA 交换（`B88` / `PbeX`）等。
#[pyclass(eq, eq_int)]
#[derive(Clone, Copy, PartialEq)]
pub enum XcComponent {
    LdaX,
    Vwn5,
    Vwn3,
    B88,
    Lyp,
    P86,
    PbeX,
    PbeC,
}

impl XcComponent {
    /// 转换成 rswfn 对应的枚举值
    fn to_core(&self) -> rswfn::rsxcs::XcComponent {
        match self {
            Self::LdaX => rswfn::rsxcs::XcComponent::LdaX,
            Self::Vwn5 => rswfn::rsxcs::XcComponent::Vwn5,
            Self::Vwn3 => rswfn::rsxcs::XcComponent::Vwn3,
            Self::B88 => rswfn::rsxcs::XcComponent::B88,
            Self::Lyp => rswfn::rsxcs::XcComponent::Lyp,
            Self::P86 => rswfn::rsxcs::XcComponent::P86,
            Self::PbeX => rswfn::rsxcs::XcComponent::PbeX,
            Self::PbeC => rswfn::rsxcs::XcComponent::PbeC,
        }
    }
}

#[pymethods]
impl XcComponent {
    /// 所有可用的泛函成分
    #[classattr]
    pub const ALL: [XcComponent; 8] = [
        XcComponent::LdaX,
        XcComponent::Vwn5,
        XcComponent::Vwn3,
        XcComponent::B88,
        XcComponent::Lyp,
        XcComponent::P86,
        XcComponent::PbeX,
        XcComponent::PbeC,
    ];

    /// 该成分是否为广义梯度近似（GGA）
    ///
    /// ## 返回
    /// GGA 成分返回 True，纯 LDA 成分返回 False
    pub fn is_gga(&self) -> bool {
        self.to_core().is_gga()
    }
}

/// 泛函的具体种类（内部使用）
#[derive(Clone, Copy, PartialEq, Eq)]
enum XcFuncKind {
    Svwn5,
    BLYP,
    BP86,
    PBE,
    PBE0,
    B3LYP,
    External(u32),
}

/// 交换关联泛函
///
/// 内置 Svwn5、BLYP、BP86、PBE、PBE0、B3LYP 六种泛函，也支持通过
/// `XcFunc.External(id)` 引用外部注册的泛函。
#[pyclass(eq)]
#[derive(Clone, Copy, PartialEq, Eq)]
pub struct XcFunc {
    kind: XcFuncKind,
}

impl XcFunc {
    /// 内部构造
    fn from_kind(kind: XcFuncKind) -> Self {
        Self { kind }
    }

    /// 转换成 rswfn 对应的枚举值
    fn to_core(&self) -> rswfn::rsxcs::XcFunc {
        match self.kind {
            XcFuncKind::Svwn5 => rswfn::rsxcs::XcFunc::Svwn5,
            XcFuncKind::BLYP => rswfn::rsxcs::XcFunc::BLYP,
            XcFuncKind::BP86 => rswfn::rsxcs::XcFunc::BP86,
            XcFuncKind::PBE => rswfn::rsxcs::XcFunc::PBE,
            XcFuncKind::PBE0 => rswfn::rsxcs::XcFunc::PBE0,
            XcFuncKind::B3LYP => rswfn::rsxcs::XcFunc::B3LYP,
            XcFuncKind::External(id) => rswfn::rsxcs::XcFunc::External(id),
        }
    }

    /// 由 rswfn 的枚举值构造
    fn from_core(core: rswfn::rsxcs::XcFunc) -> Self {
        use rswfn::rsxcs::XcFunc as Core;
        let kind = match core {
            Core::Svwn5 => XcFuncKind::Svwn5,
            Core::BLYP => XcFuncKind::BLYP,
            Core::BP86 => XcFuncKind::BP86,
            Core::PBE => XcFuncKind::PBE,
            Core::PBE0 => XcFuncKind::PBE0,
            Core::B3LYP => XcFuncKind::B3LYP,
            Core::External(id) => XcFuncKind::External(id),
        };
        Self { kind }
    }

    /// 仅根据自身保存的种类生成名称，用于 repr
    fn label(&self) -> String {
        match self.kind {
            XcFuncKind::Svwn5 => "Svwn5".to_string(),
            XcFuncKind::BLYP => "BLYP".to_string(),
            XcFuncKind::BP86 => "BP86".to_string(),
            XcFuncKind::PBE => "PBE".to_string(),
            XcFuncKind::PBE0 => "PBE0".to_string(),
            XcFuncKind::B3LYP => "B3LYP".to_string(),
            XcFuncKind::External(id) => format!("External({})", id),
        }
    }
}

#[pymethods]
impl XcFunc {
    /// 所有内置泛函
    #[classattr]
    pub const ALL: [XcFunc; 6] = [
        XcFunc {
            kind: XcFuncKind::Svwn5,
        },
        XcFunc {
            kind: XcFuncKind::BLYP,
        },
        XcFunc {
            kind: XcFuncKind::BP86,
        },
        XcFunc {
            kind: XcFuncKind::PBE,
        },
        XcFunc {
            kind: XcFuncKind::PBE0,
        },
        XcFunc {
            kind: XcFuncKind::B3LYP,
        },
    ];

    /// Svwn5 泛函（LDA）
    #[classattr]
    pub fn Svwn5() -> Self {
        Self::from_kind(XcFuncKind::Svwn5)
    }

    /// BLYP 泛函
    #[classattr]
    pub fn BLYP() -> Self {
        Self::from_kind(XcFuncKind::BLYP)
    }

    /// BP86 泛函
    #[classattr]
    pub fn BP86() -> Self {
        Self::from_kind(XcFuncKind::BP86)
    }

    /// PBE 泛函
    #[classattr]
    pub fn PBE() -> Self {
        Self::from_kind(XcFuncKind::PBE)
    }

    /// PBE0 杂化泛函
    #[classattr]
    pub fn PBE0() -> Self {
        Self::from_kind(XcFuncKind::PBE0)
    }

    /// B3LYP 杂化泛函
    #[classattr]
    pub fn B3LYP() -> Self {
        Self::from_kind(XcFuncKind::B3LYP)
    }

    /// 外部注册的泛函
    ///
    /// ## 参数
    /// - `id`: 外部泛函编号
    ///
    /// ## 返回
    /// 返回对应的 XcFunc 对象
    #[staticmethod]
    pub fn External(id: u32) -> Self {
        Self::from_kind(XcFuncKind::External(id))
    }

    /// 按名称查找泛函（不区分大小写）
    ///
    /// ## 参数
    /// - `name`: 泛函名称，如 "b3lyp"
    ///
    /// ## 返回
    /// 找到时返回 XcFunc 对象，否则返回 None
    #[staticmethod]
    pub fn from_name(name: &str) -> Option<XcFunc> {
        rswfn::rsxcs::XcFunc::from_name(name).map(Self::from_core)
    }

    /// 泛函名称
    ///
    /// ## 返回
    /// 返回泛函名称字符串
    pub fn name(&self) -> String {
        self.to_core().name().to_string()
    }

    /// 精确交换（Hartree-Fock 交换）的混合比例
    ///
    /// ## 返回
    /// 返回 0 到 1 之间的比例系数
    pub fn exact_x(&self) -> f64 {
        self.to_core().exact_x()
    }

    /// 是否为杂化泛函
    ///
    /// ## 返回
    /// 杂化泛函返回 True，纯泛函返回 False
    pub fn is_hybrid(&self) -> bool {
        self.to_core().is_hybrid()
    }

    /// 外部泛函编号
    ///
    /// ## 返回
    /// 内置泛函返回 None，外部泛函返回其编号
    pub fn external_id(&self) -> Option<u32> {
        self.to_core().external_id()
    }

    /// 返回泛函的字符串表示
    pub fn __repr__(&self) -> String {
        format!("XcFunc({})", self.label())
    }
}

/// RKS（闭壳层）泛函求值结果
#[pyclass]
pub struct RksEval {
    pub core: rswfn::rsxcs::RksEval,
}

#[pymethods]
impl RksEval {
    /// 能量密度
    #[getter(exc)]
    pub fn get_exc(&self) -> Vec<f64> {
        self.core.exc.clone()
    }

    /// 能量密度对密度 rho 的导数
    #[getter(vrho)]
    pub fn get_vrho(&self) -> Vec<f64> {
        self.core.vrho.clone()
    }

    /// 能量密度对梯度平方 sigma 的导数
    #[getter(vsigma)]
    pub fn get_vsigma(&self) -> Vec<f64> {
        self.core.vsigma.clone()
    }

    /// 返回求值结果的字符串表示
    pub fn __repr__(&self) -> String {
        format!("RksEval(npoints={})", self.core.exc.len())
    }
}

/// UKS（开壳层）泛函求值结果
#[pyclass]
pub struct UksEval {
    pub core: rswfn::rsxcs::UksEval,
}

#[pymethods]
impl UksEval {
    /// 能量密度
    #[getter(exc)]
    pub fn get_exc(&self) -> Vec<f64> {
        self.core.exc.clone()
    }

    /// 能量密度对 alpha 密度 rho_a 的导数
    #[getter(vrho_a)]
    pub fn get_vrho_a(&self) -> Vec<f64> {
        self.core.vrho_a.clone()
    }

    /// 能量密度对 beta 密度 rho_b 的导数
    #[getter(vrho_b)]
    pub fn get_vrho_b(&self) -> Vec<f64> {
        self.core.vrho_b.clone()
    }

    /// 能量密度对 sigma_aa 的导数
    #[getter(vsigma_aa)]
    pub fn get_vsigma_aa(&self) -> Vec<f64> {
        self.core.vsigma_aa.clone()
    }

    /// 能量密度对 sigma_bb 的导数
    #[getter(vsigma_bb)]
    pub fn get_vsigma_bb(&self) -> Vec<f64> {
        self.core.vsigma_bb.clone()
    }

    /// 能量密度对 sigma_ab 的导数
    #[getter(vsigma_ab)]
    pub fn get_vsigma_ab(&self) -> Vec<f64> {
        self.core.vsigma_ab.clone()
    }

    /// 返回求值结果的字符串表示
    pub fn __repr__(&self) -> String {
        format!("UksEval(npoints={})", self.core.exc.len())
    }
}

/// 交换关联泛函求值器
///
/// 封装 `rswfn::rsxcs::RsXcs`，按泛函定义把各成分的贡献加权求和。
#[pyclass]
pub struct RsXcs {
    pub core: rswfn::rsxcs::RsXcs,
}

#[pymethods]
impl RsXcs {
    /// 创建一个默认泛函（BLYP）的求值器
    ///
    /// ## 返回
    /// 返回一个 RsXcs 对象
    #[new]
    pub fn new() -> Self {
        Self {
            core: rswfn::rsxcs::RsXcs::new(),
        }
    }

    /// 按指定泛函创建求值器
    ///
    /// ## 参数
    /// - `func`: 泛函对象，如 `XcFunc.B3LYP`
    ///
    /// ## 返回
    /// 返回一个 RsXcs 对象
    #[staticmethod]
    pub fn with_func(func: PyRef<'_, XcFunc>) -> Self {
        Self {
            core: rswfn::rsxcs::RsXcs::with_func(func.to_core()),
        }
    }

    /// 当前使用的泛函
    ///
    /// ## 返回
    /// 返回 XcFunc 对象
    pub fn func(&self) -> XcFunc {
        XcFunc::from_core(self.core.func())
    }

    /// 当前泛函的精确交换比例
    ///
    /// ## 返回
    /// 返回 0 到 1 之间的比例系数
    pub fn exact_x(&self) -> f64 {
        self.core.exact_x()
    }

    /// 在闭壳层（RKS）下求值
    ///
    /// ## 参数
    /// - `rho`: 各格点的电子密度列表
    /// - `sigma`: 各格点的密度梯度平方列表
    ///
    /// ## 返回
    /// 返回 (exc, vrho, vsigma) 三元组，均为与输入等长的列表
    pub fn run(&self, rho: Vec<f64>, sigma: Vec<f64>) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        self.core.run(&rho, &sigma)
    }

    /// 在开壳层（UKS）下求值
    ///
    /// ## 参数
    /// - `rho_a`: 各格点的 alpha 电子密度列表
    /// - `rho_b`: 各格点的 beta 电子密度列表
    /// - `saa`: 各格点的 sigma_aa 列表
    /// - `sbb`: 各格点的 sigma_bb 列表
    /// - `sab`: 各格点的 sigma_ab 列表
    ///
    /// ## 返回
    /// 返回 (exc, vrho_a, vrho_b, vsigma_aa, vsigma_bb, vsigma_ab) 六元组
    pub fn run_spin(
        &self,
        rho_a: Vec<f64>,
        rho_b: Vec<f64>,
        saa: Vec<f64>,
        sbb: Vec<f64>,
        sab: Vec<f64>,
    ) -> (
        Vec<f64>,
        Vec<f64>,
        Vec<f64>,
        Vec<f64>,
        Vec<f64>,
        Vec<f64>,
    ) {
        self.core.run_spin(&rho_a, &rho_b, &saa, &sbb, &sab)
    }

    /// 返回求值器的字符串表示
    pub fn __repr__(&self) -> String {
        format!("RsXcs({})", self.func().label())
    }
}

/// 在闭壳层（RKS）下计算单个泛函成分的贡献
///
/// ## 参数
/// - `comp`: 泛函成分，如 `XcComponent.LdaX`
/// - `rho`: 各格点的电子密度列表
/// - `sigma`: 各格点的密度梯度平方列表
///
/// ## 返回
/// 返回 RksEval 对象，包含 exc、vrho、vsigma
#[pyfunction]
pub fn eval_rks(comp: PyRef<'_, XcComponent>, rho: Vec<f64>, sigma: Vec<f64>) -> RksEval {
    let core = rswfn::rsxcs::eval_rks(comp.to_core(), &rho, &sigma);
    RksEval { core }
}

/// 在开壳层（UKS）下计算单个泛函成分的贡献
///
/// ## 参数
/// - `comp`: 泛函成分，如 `XcComponent.B88`
/// - `rho_a`: 各格点的 alpha 电子密度列表
/// - `rho_b`: 各格点的 beta 电子密度列表
/// - `saa`: 各格点的 sigma_aa 列表
/// - `sbb`: 各格点的 sigma_bb 列表
/// - `sab`: 各格点的 sigma_ab 列表
///
/// ## 返回
/// 返回 UksEval 对象，包含 exc、vrho_a、vrho_b、vsigma_aa、vsigma_bb、vsigma_ab
#[pyfunction]
pub fn eval_uks(
    comp: PyRef<'_, XcComponent>,
    rho_a: Vec<f64>,
    rho_b: Vec<f64>,
    saa: Vec<f64>,
    sbb: Vec<f64>,
    sab: Vec<f64>,
) -> UksEval {
    let core = rswfn::rsxcs::eval_uks(comp.to_core(), &rho_a, &rho_b, &saa, &sbb, &sab);
    UksEval { core }
}

/// 返回所有内置泛函
///
/// ## 返回
/// 返回 XcFunc 对象列表
#[pyfunction]
pub fn supported() -> Vec<XcFunc> {
    rswfn::rsxcs::supported()
        .into_iter()
        .map(XcFunc::from_core)
        .collect()
}

/// 返回所有内置泛函的名称
///
/// ## 返回
/// 返回泛函名称列表
#[pyfunction]
pub fn supported_names() -> Vec<String> {
    rswfn::rsxcs::supported_names()
        .into_iter()
        .map(|name| name.to_string())
        .collect()
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "rsxcs")?;
    m.add_class::<XcComponent>()?;
    m.add_class::<XcFunc>()?;
    m.add_class::<RksEval>()?;
    m.add_class::<UksEval>()?;
    m.add_class::<RsXcs>()?;
    m.add_function(wrap_pyfunction!(eval_rks, &m)?)?;
    m.add_function(wrap_pyfunction!(eval_uks, &m)?)?;
    m.add_function(wrap_pyfunction!(supported, &m)?)?;
    m.add_function(wrap_pyfunction!(supported_names, &m)?)?;
    parent_module.add_submodule(&m)
}
