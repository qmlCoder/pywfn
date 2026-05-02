use pyo3::prelude::*;

mod mole;

mod atoms;
mod basis;
mod coefs;

pub use atoms::{Atom, Stm};
pub use mole::Mole;

pub use atoms::Atoms;
pub use basis::Basis;
pub use coefs::Coefs;

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "base")?;
    m.add_class::<Stm>()?;
    m.add_class::<Atom>()?;
    m.add_class::<Atoms>()?;
    m.add_class::<Basis>()?;
    m.add_class::<Coefs>()?;
    m.add_class::<Mole>()?;
    parent_module.add_submodule(&m)
}

#[pymodule]
pub mod base {
    #[pymodule_export]
    use super::Atom;
    #[pymodule_export]
    use super::Atoms;
    #[pymodule_export]
    use super::Basis;
    #[pymodule_export]
    use super::Coefs;
    #[pymodule_export]
    use super::Mole;
    #[pymodule_export]
    use super::Stm;
}
