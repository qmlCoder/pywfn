use pyo3::prelude::*;

pub mod hf;
pub mod se;
pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "qmcalc")?;
    hf::register_module(&m)?;
    se::register_module(&m)?;
    parent_module.add_submodule(&m)
}
