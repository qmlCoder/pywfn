#![allow(non_snake_case)]

use numpy::{IntoPyArray, PyArray2, PyReadonlyArray2};
use pyo3::prelude::*;
use rswfn;

/// 计算两个基函数的重叠积分
#[pyfunction]
pub fn gtf_integ(alps: [f64; 2], xyzs: [[f64; 3]; 2], lmns: [[u32; 3]; 2]) -> f64 {
    rswfn::integ::gtf_integ(&alps, &xyzs, &lmns)
}

#[pyfunction]
pub fn local_gtf_integ(
    alps: [f64; 2],
    xyzs: [[f64; 3]; 2],
    syms: [String; 2],
    stm1: PyReadonlyArray2<f64>,
    stm2: PyReadonlyArray2<f64>,
) -> f64 {
    let stm1 = stm1.as_array().to_owned();
    let stm2 = stm2.as_array().to_owned();
    let syms: [&str; 2] = [syms[0].as_str(), syms[1].as_str()];
    rswfn::integ::local_gtf_integ(&alps, &xyzs, &syms, &[stm1, stm2])
}

#[pyfunction]
pub fn calc_smat(
    py: Python,
    alps: Vec<Vec<f64>>,
    coes: Vec<Vec<f64>>,
    lmns: Vec<[u32; 3]>,
    xyzs: Vec<[f64; 3]>,
) -> Py<PyArray2<f64>> {
    let smat = rswfn::integ::calc_smat(&alps, &coes, &lmns, &xyzs);
    smat.into_pyarray(py).unbind()
}

#[pyfunction]
pub fn rys_roots(n: usize, alp: f64) -> Vec<[f64; 2]> {
    rswfn::integ::ry::rys_roots(n, alp)
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "integ")?; // 为当前的rs文件创建一个子模块
    m.add_function(wrap_pyfunction!(gtf_integ, &m)?)?;
    m.add_function(wrap_pyfunction!(calc_smat, &m)?)?;
    m.add_function(wrap_pyfunction!(local_gtf_integ, &m)?)?;
    m.add_function(wrap_pyfunction!(rys_roots, &m)?)?;
    parent_module.add_submodule(&m) // 将子模块添加到父模块中
}
