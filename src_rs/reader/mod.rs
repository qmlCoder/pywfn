mod cub;
mod fch;
mod gjf;
mod log;
mod mde;
mod mol;
mod non;
mod sdf;
mod xyz;

use pyo3::prelude::*;

pub use cub::CubReader;
pub use fch::FchReader;
pub use gjf::GjfReader;
pub use log::LogReader;
pub use mde::MdeReader;
pub use mol::MolReader;
pub use non::NonReader;
pub use sdf::SdfReader;
pub use xyz::XyzReader;

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "reader")?;
    m.add_class::<FchReader>()?;
    m.add_class::<GjfReader>()?;
    m.add_class::<LogReader>()?;
    m.add_class::<MolReader>()?;
    m.add_class::<NonReader>()?;
    m.add_class::<SdfReader>()?;
    m.add_class::<XyzReader>()?;
    m.add_class::<MdeReader>()?;
    m.add_class::<CubReader>()?;
    parent_module.add_submodule(&m)
}

#[pymodule]
pub mod reader {
    #[pymodule_export]
    use super::fch::FchReader;
    #[pymodule_export]
    use super::gjf::GjfReader;
    #[pymodule_export]
    use super::log::LogReader;
    #[pymodule_export]
    use super::mde::MdeReader;
    #[pymodule_export]
    use super::mol::MolReader;
    #[pymodule_export]
    use super::non::NonReader;
    #[pymodule_export]
    use super::sdf::SdfReader;
    #[pymodule_export]
    use super::xyz::XyzReader;
}
