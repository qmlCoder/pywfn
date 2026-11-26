use pyo3::prelude::*;

mod mole;

mod atoms;
mod basis;
mod coefs;
mod vect;

pub use mole::Mole;

pub use atoms::{Atom, Atoms, Stm};
pub use basis::{Atos, Basis, Shl, Shls};
pub use coefs::Coefs;
pub use vect::Vec3;

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "base")?;
    m.add_class::<Atoms>()?;
    m.add_class::<Stm>()?;
    m.add_class::<Atom>()?;
    m.add_class::<Shl>()?;
    m.add_class::<Shls>()?;

    m.add_class::<Basis>()?;
    m.add_class::<Atos>()?;

    m.add_class::<Coefs>()?;

    m.add_class::<Mole>()?;
    m.add_class::<Vec3>()?;
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
    #[pymodule_export]
    use super::Vec3;
}
