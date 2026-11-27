use crate::base::Mole;
use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use rswfn;

/// π/σ/πσ 能量分解结果（闭壳层，杂化泛函）。
///
/// 把完整电子能量按密度拆分 `D = D_π + D_σ` 做加法分解，满足
/// `E_π + E_σ + E_πσ = E_elec`。
#[pyclass]
#[derive(Clone, Copy)]
pub struct PiEnergy {
    core: rswfn::props::moleprop::energy::PiEnergy,
}

#[pymethods]
impl PiEnergy {
    /// π 自能量 `E_π`
    #[getter]
    pub fn e_pi(&self) -> f64 {
        self.core.e_pi
    }
    /// σ 自能量 `E_σ`
    #[getter]
    pub fn e_sigma(&self) -> f64 {
        self.core.e_sigma
    }
    /// π-σ 交叉项 `E_πσ`（含 XC 交叉/残差）
    #[getter]
    pub fn e_cross(&self) -> f64 {
        self.core.e_cross
    }
    /// 交叉项的库仑部分 `Tr(D_π J_σ)`
    #[getter]
    pub fn e_cross_coul(&self) -> f64 {
        self.core.e_cross_coul
    }
    /// 交叉项的精确交换部分 `−aₓ/2·Tr(D_π K_σ)`
    #[getter]
    pub fn e_cross_exch(&self) -> f64 {
        self.core.e_cross_exch
    }
    /// 交叉项的 XC 残差
    #[getter]
    pub fn e_cross_xc(&self) -> f64 {
        self.core.e_cross_xc
    }
    /// 完整电子能量 `E_elec`
    #[getter]
    pub fn e_elec(&self) -> f64 {
        self.core.e_elec
    }
    /// 核-核排斥能 `E_nuc`
    #[getter]
    pub fn e_nuc(&self) -> f64 {
        self.core.e_nuc
    }
    /// 总能量 `E_total = E_elec + E_nuc`
    #[getter]
    pub fn e_total(&self) -> f64 {
        self.core.e_total
    }

    /// π 对电子能量的可加贡献 `E_π + E_πσ`
    pub fn pi(&self) -> f64 {
        self.core.pi()
    }
}

/// 能量计算器。
///
/// 基于分子轨道信息，对完整电子能量做 π/σ/πσ 分解。
#[pyclass]
pub struct Calculator {
    pub mole: Mole,
}

impl Calculator {
    // 提取公共的 calculator 创建逻辑
    fn caler(&self) -> rswfn::props::moleprop::energy::Calculator<'_> {
        rswfn::props::moleprop::energy::Calculator::new(&self.mole.core)
    }
}

#[pymethods]
impl Calculator {
    /// 创建一个新的能量计算器。
    ///
    /// # 参数
    /// * `mole` - 分子对象（需来自高斯 log/out，且为闭壳层、含轨道系数）。
    #[new]
    pub fn new(mole: Mole) -> Self {
        Self { mole }
    }

    /// POCV 方法的 π/σ/πσ 能量分解。
    pub fn pi_pocv_parts(&self) -> PyResult<PiEnergy> {
        self.caler()
            .pi_pocv_parts()
            .map(|core| PiEnergy { core })
            .map_err(|e| PyRuntimeError::new_err(e))
    }

    /// MOCV 方法的 π/σ/πσ 能量分解。
    pub fn pi_mocv_parts(&self) -> PyResult<PiEnergy> {
        self.caler()
            .pi_mocv_parts()
            .map(|core| PiEnergy { core })
            .map_err(|e| PyRuntimeError::new_err(e))
    }

    /// 使用 POCV 方法计算的 π 对电子能量的贡献 `E_π + E_πσ`。
    pub fn pi_pocv(&self) -> PyResult<f64> {
        self.caler()
            .pi_pocv()
            .map_err(|e| PyRuntimeError::new_err(e))
    }

    /// 使用 MOCV 方法计算的 π 对电子能量的贡献 `E_π + E_πσ`。
    pub fn pi_mocv(&self) -> PyResult<f64> {
        self.caler()
            .pi_mocv()
            .map_err(|e| PyRuntimeError::new_err(e))
    }
}

/// 在父模块中注册 `energy` 子模块。
pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "energy")?;
    m.add_class::<Calculator>()?;
    m.add_class::<PiEnergy>()?;
    parent_module.add_submodule(&m)
}
