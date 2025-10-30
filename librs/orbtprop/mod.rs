use pyo3::prelude::*;

pub mod obtmat;
pub mod obtval;

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "orbtprop")?;
    obtmat::register_module(&m)?;
    obtval::register_module(&m)?;
    parent_module.add_submodule(&m)
}
