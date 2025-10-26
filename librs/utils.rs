#![allow(non_snake_case)]

use ndarray::Array2;
use pyo3::prelude::*;
use rswfn;

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "utils")?;
    parent_module.add_submodule(&m)
}
