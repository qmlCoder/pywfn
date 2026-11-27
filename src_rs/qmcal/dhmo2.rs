use crate::base::{Mole, Stm};
use pyo3::prelude::*;
use rswfn;

/// 计算两个位于 x 轴上的基函数的重叠积分（距离 > 3.0 时返回 0）
#[pyfunction]
pub fn get_integ_x(zi: usize, zj: usize, pi: usize, pj: usize, dist: f64) -> f64 {
    rswfn::qmcal::dhmo2::get_integ_x(zi, zj, pi, pj, dist)
}

/// 计算两个任意朝向的 Pz 轨道的重叠积分
#[pyfunction]
pub fn get_integ_d(
    zi: usize,
    zj: usize,
    di: [f64; 3],
    dj: [f64; 3],
    stm: Stm,
    dist: f64,
) -> f64 {
    let di = rswfn::base::Vec3::from_arr(di);
    let dj = rswfn::base::Vec3::from_arr(dj);
    rswfn::qmcal::dhmo2::get_integ_d(zi, zj, &di, &dj, stm.core, dist)
}

/// dHMO 方法版本 2，计算分子中的所有原子轨道（输出到日志）
#[pyfunction]
pub fn dhmo_calc_v2(mole: Mole) {
    rswfn::qmcal::dhmo2::dhmo_calc_v2(&mole.core)
}

/// 在父模块中注册 `dhmo2` 子模块。
pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "dhmo2")?;
    m.add_function(wrap_pyfunction!(get_integ_x, &m)?)?;
    m.add_function(wrap_pyfunction!(get_integ_d, &m)?)?;
    m.add_function(wrap_pyfunction!(dhmo_calc_v2, &m)?)?;
    parent_module.add_submodule(&m)
}
