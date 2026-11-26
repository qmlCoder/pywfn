use pyo3::prelude::*;

mod cub;
mod dmo;
mod fch;
mod gjf;
mod obj;
mod xyz;

pub use cub::CubWriter;
pub use dmo::DmoWriter;
pub use fch::FchWriter;
pub use gjf::GjfWriter;
pub use obj::ObjWriter;
pub use xyz::XyzWriter;

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "writer")?;
    m.add_class::<CubWriter>()?;
    m.add_class::<DmoWriter>()?;
    m.add_class::<FchWriter>()?;
    m.add_class::<GjfWriter>()?;
    m.add_class::<ObjWriter>()?;
    m.add_class::<XyzWriter>()?;
    parent_module.add_submodule(&m)
}

#[pymodule]
pub mod writer {
    #[pymodule_export]
    use super::CubWriter;
    #[pymodule_export]
    use super::DmoWriter;
    #[pymodule_export]
    use super::FchWriter;
    #[pymodule_export]
    use super::GjfWriter;
    #[pymodule_export]
    use super::ObjWriter;
    #[pymodule_export]
    use super::XyzWriter;
}
