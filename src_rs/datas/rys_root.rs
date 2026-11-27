use pyo3::prelude::*;
use rswfn;

/// 把固定长度的一维常量数组转换为列表
fn poly(vals: &[f64; 496]) -> Vec<f64> {
    vals.to_vec()
}

/// 大 X 区间的根
#[pyfunction]
fn POLY_LARGEX_RT() -> Vec<f64> {
    poly(&rswfn::datas::rys_root::POLY_LARGEX_RT)
}

/// 大 X 区间的权重
#[pyfunction]
fn POLY_LARGEX_WW() -> Vec<f64> {
    poly(&rswfn::datas::rys_root::POLY_LARGEX_WW)
}

/// 小 X 区间的第一组根
#[pyfunction]
fn POLY_SMALLX_R0() -> Vec<f64> {
    poly(&rswfn::datas::rys_root::POLY_SMALLX_R0)
}

/// 小 X 区间的第二组根
#[pyfunction]
fn POLY_SMALLX_R1() -> Vec<f64> {
    poly(&rswfn::datas::rys_root::POLY_SMALLX_R1)
}

/// 小 X 区间的第一组权重
#[pyfunction]
fn POLY_SMALLX_W0() -> Vec<f64> {
    poly(&rswfn::datas::rys_root::POLY_SMALLX_W0)
}

/// 小 X 区间的第二组权重
#[pyfunction]
fn POLY_SMALLX_W1() -> Vec<f64> {
    poly(&rswfn::datas::rys_root::POLY_SMALLX_W1)
}

/// 计算 Rys 1 阶多项式的根和权重
///
/// ## 参数
/// - `X`: Rys 多项式参数
///
/// ## 返回
/// 返回 (roots, weights) 元组，分别为根和权重列表
#[pyfunction]
pub fn rys_root1(X: f64) -> (Vec<f64>, Vec<f64>) {
    rswfn::datas::rys_root::rys_root1(X)
}

/// 计算 Rys 2 阶多项式的根和权重
///
/// ## 参数
/// - `X`: Rys 多项式参数
///
/// ## 返回
/// 返回 (roots, weights) 元组，分别为根和权重列表
#[pyfunction]
pub fn rys_root2(X: f64) -> (Vec<f64>, Vec<f64>) {
    rswfn::datas::rys_root::rys_root2(X)
}

/// 计算 Rys 3 阶多项式的根和权重
///
/// ## 参数
/// - `X`: Rys 多项式参数
///
/// ## 返回
/// 返回 (roots, weights) 元组，分别为根和权重列表
#[pyfunction]
pub fn rys_root3(X: f64) -> (Vec<f64>, Vec<f64>) {
    rswfn::datas::rys_root::rys_root3(X)
}

/// 计算 Rys 4 阶多项式的根和权重
///
/// ## 参数
/// - `X`: Rys 多项式参数
///
/// ## 返回
/// 返回 (roots, weights) 元组，分别为根和权重列表
#[pyfunction]
pub fn rys_root4(X: f64) -> (Vec<f64>, Vec<f64>) {
    rswfn::datas::rys_root::rys_root4(X)
}

/// 计算 Rys 5 阶多项式的根和权重
///
/// ## 参数
/// - `X`: Rys 多项式参数
///
/// ## 返回
/// 返回 (roots, weights) 元组，分别为根和权重列表
#[pyfunction]
pub fn rys_root5(X: f64) -> (Vec<f64>, Vec<f64>) {
    rswfn::datas::rys_root::rys_root5(X)
}

/// 计算任意阶 Rys 多项式的根和权重
///
/// ## 参数
/// - `n`: 阶数
/// - `x`: Rys 多项式参数
///
/// ## 返回
/// 返回每点的 (根, 权重) 列表
#[pyfunction]
pub fn rys_roots(n: usize, x: f64) -> Vec<[f64; 2]> {
    rswfn::datas::rys_root::rys_roots(n, x)
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "rys_root")?;
    m.add_function(wrap_pyfunction!(POLY_LARGEX_RT, &m)?)?;
    m.add_function(wrap_pyfunction!(POLY_LARGEX_WW, &m)?)?;
    m.add_function(wrap_pyfunction!(POLY_SMALLX_R0, &m)?)?;
    m.add_function(wrap_pyfunction!(POLY_SMALLX_R1, &m)?)?;
    m.add_function(wrap_pyfunction!(POLY_SMALLX_W0, &m)?)?;
    m.add_function(wrap_pyfunction!(POLY_SMALLX_W1, &m)?)?;
    m.add_function(wrap_pyfunction!(rys_root1, &m)?)?;
    m.add_function(wrap_pyfunction!(rys_root2, &m)?)?;
    m.add_function(wrap_pyfunction!(rys_root3, &m)?)?;
    m.add_function(wrap_pyfunction!(rys_root4, &m)?)?;
    m.add_function(wrap_pyfunction!(rys_root5, &m)?)?;
    m.add_function(wrap_pyfunction!(rys_roots, &m)?)?;
    parent_module.add_submodule(&m)
}
