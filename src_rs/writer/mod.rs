use pyo3::prelude::*;

mod cub;
mod gjf;

pub use cub::CubWriter;
pub use gjf::GjfWriter;

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "writer")?;
    m.add_class::<GjfWriter>()?;
    m.add_class::<CubWriter>()?;
    parent_module.add_submodule(&m)
}

#[pymodule]
pub mod writer {
    #[pymodule_export]
    use super::GjfWriter;
}
