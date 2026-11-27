use pyo3::prelude::*;
use rswfn;

/// 获取 HMO（Hückel 分子轨道）积分参数
///
/// ## 参数
/// - `Zi`: 原子 i 的电荷数
/// - `Zj`: 原子 j 的电荷数
/// - `Pi`: 原子 i 的主量子数
/// - `Pj`: 原子 j 的主量子数
///
/// ## 返回
/// 返回 59 组参数，每组 4 个值
#[pyfunction]
fn get_paras(Zi: usize, Zj: usize, Pi: usize, Pj: usize) -> Vec<[f64; 4]> {
    rswfn::datas::dhmo2::get_paras(Zi, Zj, Pi, Pj).to_vec()
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "dhmo2")?;
    m.add_function(wrap_pyfunction!(get_paras, &m)?)?;
    parent_module.add_submodule(&m)
}
