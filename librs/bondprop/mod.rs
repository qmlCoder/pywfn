use pyo3::prelude::*;

pub mod direction;
pub mod order;

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "bondprop")?;
    direction::register_module(&m)?;
    order::register_module(&m)?;
    parent_module.add_submodule(&m)
}
