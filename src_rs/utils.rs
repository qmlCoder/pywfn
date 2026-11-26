#![allow(non_snake_case)]

use ndarray::Array2;
use pyo3::prelude::*;
use rswfn;

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "utils")?;
    m.add_function(wrap_pyfunction!(rswfn_version, &m)?)?;
    parent_module.add_submodule(&m)
}

/// 返回底层依赖的 rswfn 库版本号。
///
/// 该值在构建时从 Cargo.lock 中读取，反映实际编译链接的 rswfn 版本。
#[pyfunction]
pub fn rswfn_version() -> &'static str {
    env!("RSWFN_VERSION")
}

pub fn arr2vec(arr: &Array2<f64>) -> Vec<Vec<f64>> {
    let nrow = arr.nrows();
    let ncol = arr.ncols();
    let mut vec = vec![vec![0.0; ncol]; nrow];
    for i in 0..nrow {
        for j in 0..ncol {
            vec[i][j] = arr[[i, j]];
        }
    }
    vec
}
