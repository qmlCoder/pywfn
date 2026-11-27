use pyo3::prelude::*;
use rswfn;

pub mod dhmo;
pub mod dhmo2;
pub mod diis;
pub mod guess;
pub mod hf;
pub mod ks;
pub mod lks;
pub mod se;

/// Slater 型原子轨道的指数换算函数
///
/// 按 Slater 规则把标度参数换算成对应的轨道指数，供 dHMO（双中心 Hückel）
/// 计算构造 Slater 型基函数时使用。
///
/// ## 参数
/// - `rho`: 标度参数（无量纲）
///
/// ## 返回
/// 换算后的 Slater 指数
#[pyfunction]
#[allow(non_snake_case)]
pub fn SlaterX(rho: f64) -> f64 {
    rswfn::qmcal::SlaterX(rho)
}

/// 在父模块中注册 `qmcal` 子模块。
pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "qmcal")?;
    m.add_function(wrap_pyfunction!(SlaterX, &m)?)?;
    hf::register_module(&m)?;
    se::register_module(&m)?;
    ks::register_module(&m)?;
    lks::register_module(&m)?;
    diis::register_module(&m)?;
    dhmo::register_module(&m)?;
    dhmo2::register_module(&m)?;
    guess::register_module(&m)?;
    parent_module.add_submodule(&m)
}
