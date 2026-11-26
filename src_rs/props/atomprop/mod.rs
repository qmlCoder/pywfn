use pyo3::prelude::*;

pub mod activity;
pub mod charge;
pub mod direction;

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "atomprop")?;
    charge::register_module(&m)?;
    activity::register_module(&m)?;
    direction::register_module(&m)?;
    parent_module.add_submodule(&m)
}
