#![allow(non_snake_case)]

use pyo3::prelude::*;

mod utils;

pub mod base;
pub mod datas;
pub mod integ;
pub mod march;
pub mod maths;
pub mod matrix;
pub mod qmcalc;
pub mod reader;
pub mod space;
pub mod tools;
pub mod writer;

pub mod atomprop;
pub mod bondprop;
pub mod fragprop;
pub mod gridprop;
pub mod moleprop;
pub mod orbtprop;

#[pymodule]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    base::register_module(m)?;
    space::register_module(m)?;
    reader::register_module(m)?;
    writer::register_module(m)?;
    utils::register_module(m)?;
    march::register_module(m)?;
    maths::register_module(m)?;
    integ::register_module(m)?;
    matrix::register_module(m)?;
    tools::register_module(m)?;
    datas::register_module(m)?;

    atomprop::register_module(m)?;
    bondprop::register_module(m)?;
    fragprop::register_module(m)?;
    gridprop::register_module(m)?;
    moleprop::register_module(m)?;
    orbtprop::register_module(m)?;
    Ok(())
}

// #[pymodule]
// mod core {

//     use pyo3::pymodule;

//     #[pymodule_export]
//     use super::base::base;

//     #[pymodule_export]
//     use super::reader::reader;
//     #[pymodule_export]
//     use super::writer::writer;

//     #[pymodule_export]
//     use super::atomprop::atomprop;
//     #[pymodule_export]
//     use super::bondprop::bondprop;
//     #[pymodule_export]
//     use super::fragprop::fragprop;
//     #[pymodule_export]
//     use super::gridprop::gridprop;
//     #[pymodule_export]
//     use super::moleprop::moleprop;
//     #[pymodule_export]
//     use super::orbtprop::orbtprop;
// }
