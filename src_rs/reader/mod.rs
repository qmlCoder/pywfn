mod cub;
mod dmo;
mod fch;
mod gjf;
mod inp;
mod log;
mod mde;
mod mol;
mod non;
mod out;
mod sdf;
mod xyz;

use pyo3::prelude::*;

pub use cub::CubReader;
pub use dmo::DmoReader;
pub use fch::FchReader;
pub use gjf::GjfReader;
pub use inp::InpReader;
pub use log::LogReader;
pub use mde::MdeReader;
pub use mol::MolReader;
pub use non::NonReader;
pub use out::OutReader;
pub use sdf::SdfReader;
pub use xyz::XyzReader;

use rswfn;

#[pyclass]
pub struct Nics {
    core: rswfn::reader::Nics,
}

#[pymethods]
impl Nics {
    #[getter]
    fn iso(&self) -> f64 {
        self.core.iso
    }

    #[getter]
    fn ani(&self) -> f64 {
        self.core.ani
    }

    #[getter]
    fn mat(&self) -> [f64; 9] {
        self.core.mat
    }

    #[getter]
    fn eig(&self) -> [f64; 3] {
        self.core.eig
    }
}

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
    m.add_class::<DmoReader>()?;
    m.add_class::<InpReader>()?;
    m.add_class::<OutReader>()?;
    m.add_class::<Nics>()?;
    parent_module.add_submodule(&m)
}
