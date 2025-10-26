use pyo3::prelude::*;

pub mod aromacity;
pub mod population;

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "moleprop")?;
    population::register_module(&m)?;
    aromacity::register_module(&m)?;
    parent_module.add_submodule(&m)
}
