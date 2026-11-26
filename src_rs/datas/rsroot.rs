use pyo3::prelude::*;
use rswfn;

/// 计算 Rys 1 阶多项式的根和权重
///
/// ## 参数
/// - `X`: Rys 多项式参数
///
/// ## 返回
/// 返回 (roots, weights) 元组，分别为根和权重列表
#[pyfunction]
pub fn rys_root1(X: f64) -> PyResult<(Vec<f64>, Vec<f64>)> {
    Ok(rswfn::datas::rys_root::rys_root1(X))
}

/// 计算 Rys 2 阶多项式的根和权重
///
/// ## 参数
/// - `X`: Rys 多项式参数
///
/// ## 返回
/// 返回 (roots, weights) 元组，分别为根和权重列表
#[pyfunction]
pub fn rys_root2(X: f64) -> PyResult<(Vec<f64>, Vec<f64>)> {
    Ok(rswfn::datas::rys_root::rys_root2(X))
}

/// 计算 Rys 3 阶多项式的根和权重
///
/// ## 参数
/// - `X`: Rys 多项式参数
///
/// ## 返回
/// 返回 (roots, weights) 元组，分别为根和权重列表
#[pyfunction]
pub fn rys_root3(X: f64) -> PyResult<(Vec<f64>, Vec<f64>)> {
    Ok(rswfn::datas::rys_root::rys_root3(X))
}

/// 计算 Rys 4 阶多项式的根和权重
///
/// ## 参数
/// - `X`: Rys 多项式参数
///
/// ## 返回
/// 返回 (roots, weights) 元组，分别为根和权重列表
#[pyfunction]
pub fn rys_root4(X: f64) -> PyResult<(Vec<f64>, Vec<f64>)> {
    Ok(rswfn::datas::rys_root::rys_root4(X))
}

/// 计算 Rys 5 阶多项式的根和权重
///
/// ## 参数
/// - `X`: Rys 多项式参数
///
/// ## 返回
/// 返回 (roots, weights) 元组，分别为根和权重列表
#[pyfunction]
pub fn rys_root5(X: f64) -> PyResult<(Vec<f64>, Vec<f64>)> {
    Ok(rswfn::datas::rys_root::rys_root5(X))
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "rsroot")?;
    m.add_function(wrap_pyfunction!(rys_root1, &m)?)?;
    m.add_function(wrap_pyfunction!(rys_root2, &m)?)?;
    m.add_function(wrap_pyfunction!(rys_root3, &m)?)?;
    m.add_function(wrap_pyfunction!(rys_root4, &m)?)?;
    m.add_function(wrap_pyfunction!(rys_root5, &m)?)?;
    parent_module.add_submodule(&m)
}
