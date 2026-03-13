use pyo3::prelude::*;

mod gjf;

pub use gjf::GjfWriter;

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "writer")?;
    m.add_class::<GjfWriter>()?;
    parent_module.add_submodule(&m)
}
