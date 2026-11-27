use pyo3::prelude::*;
use rswfn;

/// 把固定长度的二维常量数组转换为嵌套列表
fn arrs<const N: usize>(vals: &[[f64; 50]; N]) -> Vec<Vec<f64>> {
    vals.iter().map(|row| row.to_vec()).collect()
}

/// Rys 多项式求积的 alpha 参数
#[pyfunction]
fn RYS_ALPS() -> Vec<f64> {
    rswfn::datas::roots::RYS_ALPS.to_vec()
}

/// 1 点 Rys 求积的根
#[pyfunction]
fn RYS_ROOTS_1() -> Vec<Vec<f64>> {
    arrs(&rswfn::datas::roots::RYS_ROOTS_1)
}

/// 2 点 Rys 求积的根
#[pyfunction]
fn RYS_ROOTS_2() -> Vec<Vec<f64>> {
    arrs(&rswfn::datas::roots::RYS_ROOTS_2)
}

/// 3 点 Rys 求积的根
#[pyfunction]
fn RYS_ROOTS_3() -> Vec<Vec<f64>> {
    arrs(&rswfn::datas::roots::RYS_ROOTS_3)
}

/// 4 点 Rys 求积的根
#[pyfunction]
fn RYS_ROOTS_4() -> Vec<Vec<f64>> {
    arrs(&rswfn::datas::roots::RYS_ROOTS_4)
}

/// 5 点 Rys 求积的根
#[pyfunction]
fn RYS_ROOTS_5() -> Vec<Vec<f64>> {
    arrs(&rswfn::datas::roots::RYS_ROOTS_5)
}

/// 6 点 Rys 求积的根
#[pyfunction]
fn RYS_ROOTS_6() -> Vec<Vec<f64>> {
    arrs(&rswfn::datas::roots::RYS_ROOTS_6)
}

/// 7 点 Rys 求积的根
#[pyfunction]
fn RYS_ROOTS_7() -> Vec<Vec<f64>> {
    arrs(&rswfn::datas::roots::RYS_ROOTS_7)
}

/// 8 点 Rys 求积的根
#[pyfunction]
fn RYS_ROOTS_8() -> Vec<Vec<f64>> {
    arrs(&rswfn::datas::roots::RYS_ROOTS_8)
}

/// 9 点 Rys 求积的根
#[pyfunction]
fn RYS_ROOTS_9() -> Vec<Vec<f64>> {
    arrs(&rswfn::datas::roots::RYS_ROOTS_9)
}

/// 1 点 Rys 求积的权重
#[pyfunction]
fn RYS_WEITS_1() -> Vec<Vec<f64>> {
    arrs(&rswfn::datas::roots::RYS_WEITS_1)
}

/// 2 点 Rys 求积的权重
#[pyfunction]
fn RYS_WEITS_2() -> Vec<Vec<f64>> {
    arrs(&rswfn::datas::roots::RYS_WEITS_2)
}

/// 3 点 Rys 求积的权重
#[pyfunction]
fn RYS_WEITS_3() -> Vec<Vec<f64>> {
    arrs(&rswfn::datas::roots::RYS_WEITS_3)
}

/// 4 点 Rys 求积的权重
#[pyfunction]
fn RYS_WEITS_4() -> Vec<Vec<f64>> {
    arrs(&rswfn::datas::roots::RYS_WEITS_4)
}

/// 5 点 Rys 求积的权重
#[pyfunction]
fn RYS_WEITS_5() -> Vec<Vec<f64>> {
    arrs(&rswfn::datas::roots::RYS_WEITS_5)
}

/// 6 点 Rys 求积的权重
#[pyfunction]
fn RYS_WEITS_6() -> Vec<Vec<f64>> {
    arrs(&rswfn::datas::roots::RYS_WEITS_6)
}

/// 7 点 Rys 求积的权重
#[pyfunction]
fn RYS_WEITS_7() -> Vec<Vec<f64>> {
    arrs(&rswfn::datas::roots::RYS_WEITS_7)
}

/// 8 点 Rys 求积的权重
#[pyfunction]
fn RYS_WEITS_8() -> Vec<Vec<f64>> {
    arrs(&rswfn::datas::roots::RYS_WEITS_8)
}

/// 9 点 Rys 求积的权重
#[pyfunction]
fn RYS_WEITS_9() -> Vec<Vec<f64>> {
    arrs(&rswfn::datas::roots::RYS_WEITS_9)
}

/// 打印 Rys 求积数据表
#[pyfunction]
fn dub_roots() {
    rswfn::datas::roots::dub_roots()
}

/// 计算 Hermite 求积的根和权重
///
/// ## 参数
/// - `n`: 求积点数
///
/// ## 返回
/// 返回每点的 (根, 权重) 列表
#[pyfunction]
fn hermite(n: usize) -> Vec<[f64; 2]> {
    rswfn::datas::roots::hermite(n)
}

/// 计算 Legendre 求积的根和权重
///
/// ## 参数
/// - `n`: 求积点数
///
/// ## 返回
/// 返回每点的 (根, 权重) 列表
#[pyfunction]
fn legendre(n: usize) -> Vec<[f64; 2]> {
    rswfn::datas::roots::legendre(n)
}

/// 计算 Rys 求积的根和权重
///
/// ## 参数
/// - `n`: 求积点数
/// - `alp`: Rys 参数
///
/// ## 返回
/// 返回每点的 (根, 权重) 列表
#[pyfunction]
fn rys_roots(n: usize, alp: f64) -> Vec<[f64; 2]> {
    rswfn::datas::roots::rys_roots(n, alp)
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "roots")?;
    m.add_function(wrap_pyfunction!(RYS_ALPS, &m)?)?;
    m.add_function(wrap_pyfunction!(RYS_ROOTS_1, &m)?)?;
    m.add_function(wrap_pyfunction!(RYS_ROOTS_2, &m)?)?;
    m.add_function(wrap_pyfunction!(RYS_ROOTS_3, &m)?)?;
    m.add_function(wrap_pyfunction!(RYS_ROOTS_4, &m)?)?;
    m.add_function(wrap_pyfunction!(RYS_ROOTS_5, &m)?)?;
    m.add_function(wrap_pyfunction!(RYS_ROOTS_6, &m)?)?;
    m.add_function(wrap_pyfunction!(RYS_ROOTS_7, &m)?)?;
    m.add_function(wrap_pyfunction!(RYS_ROOTS_8, &m)?)?;
    m.add_function(wrap_pyfunction!(RYS_ROOTS_9, &m)?)?;
    m.add_function(wrap_pyfunction!(RYS_WEITS_1, &m)?)?;
    m.add_function(wrap_pyfunction!(RYS_WEITS_2, &m)?)?;
    m.add_function(wrap_pyfunction!(RYS_WEITS_3, &m)?)?;
    m.add_function(wrap_pyfunction!(RYS_WEITS_4, &m)?)?;
    m.add_function(wrap_pyfunction!(RYS_WEITS_5, &m)?)?;
    m.add_function(wrap_pyfunction!(RYS_WEITS_6, &m)?)?;
    m.add_function(wrap_pyfunction!(RYS_WEITS_7, &m)?)?;
    m.add_function(wrap_pyfunction!(RYS_WEITS_8, &m)?)?;
    m.add_function(wrap_pyfunction!(RYS_WEITS_9, &m)?)?;
    m.add_function(wrap_pyfunction!(dub_roots, &m)?)?;
    m.add_function(wrap_pyfunction!(hermite, &m)?)?;
    m.add_function(wrap_pyfunction!(legendre, &m)?)?;
    m.add_function(wrap_pyfunction!(rys_roots, &m)?)?;
    parent_module.add_submodule(&m)
}
