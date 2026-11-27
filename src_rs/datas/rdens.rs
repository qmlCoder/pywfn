use std::collections::HashMap;

use pyo3::prelude::*;
use rswfn;

/// 获取指定原子的径向电子密度（方法一）
///
/// ## 参数
/// - `atomic`: 原子序数
/// - `dist`: 距离原子核的距离
#[pyfunction]
fn get_rad_dens_v1(atomic: usize, dist: f64) -> f64 {
    rswfn::datas::rdens::get_rad_dens_v1(atomic, dist)
}

/// 获取指定原子的径向电子密度（方法二）
///
/// ## 参数
/// - `atomic`: 原子序数
/// - `dist`: 距离原子核的距离
#[pyfunction]
fn get_rad_dens_v2(atomic: usize, dist: f64) -> f64 {
    rswfn::datas::rdens::get_rad_dens_v2(atomic, dist)
}

/// 原子密度的系数表
#[pyfunction]
fn ADP() -> Vec<Vec<f64>> {
    rswfn::datas::rdens::ADP.clone()
}

/// 各元素的径向密度表
#[pyfunction]
fn RAD_DENS() -> HashMap<usize, Vec<f64>> {
    rswfn::datas::rdens::RAD_DENS.clone()
}

/// 径向密度的距离网格
#[pyfunction]
fn RAD_DIST() -> Vec<f64> {
    rswfn::datas::rdens::RAD_DIST.clone()
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "rdens")?;
    m.add_function(wrap_pyfunction!(ADP, &m)?)?;
    m.add_function(wrap_pyfunction!(RAD_DENS, &m)?)?;
    m.add_function(wrap_pyfunction!(RAD_DIST, &m)?)?;
    m.add_function(wrap_pyfunction!(get_rad_dens_v1, &m)?)?;
    m.add_function(wrap_pyfunction!(get_rad_dens_v2, &m)?)?;
    parent_module.add_submodule(&m)
}
