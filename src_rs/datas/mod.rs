mod btrans;
mod rsroot;
use pyo3::prelude::*;

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "datas")?;
    btrans::register_module(&m)?;
    rsroot::register_module(&m)?;
    parent_module.add_submodule(&m)
}
