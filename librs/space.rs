#![allow(non_snake_case)]

use pyo3::prelude::*;
use rayon::prelude::*;
use rswfn;

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "space")?;
    parent_module.add_submodule(&m)
}
