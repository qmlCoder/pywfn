use pyo3::prelude::*;

mod atom;
mod mole;

mod basis;
mod coefs;
mod geome;

pub use atom::Atom;
pub use mole::Mole;

pub use basis::Basis;
pub use coefs::Coefs;
pub use geome::Geome;

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "base")?;
    m.add_class::<Mole>()?;
    m.add_class::<Atom>()?;
    m.add_class::<Geome>()?;
    m.add_class::<Basis>()?;
    m.add_class::<Coefs>()?;
    parent_module.add_submodule(&m)
}
