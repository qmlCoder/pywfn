#![allow(non_snake_case)]

use ndarray::Array2;
use pyo3::prelude::*;
use rswfn;

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "utils")?;
    parent_module.add_submodule(&m)
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
