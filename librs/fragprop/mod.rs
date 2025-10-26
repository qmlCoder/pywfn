use pyo3::prelude::*;

pub mod energy;

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "fragprop")?;
    energy::register_module(&m)?;
    energy::register_module(&m)?;
    parent_module.add_submodule(&m)
}
