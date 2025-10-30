#![allow(non_snake_case)]

use pyo3::prelude::*;

pub mod base;
pub mod integ;
pub mod march;
pub mod maths;
pub mod matrix;
pub mod reader;
pub mod space;
pub mod tools;
pub mod utils;

pub mod atomprop;
pub mod bondprop;
pub mod fragprop;
pub mod gridprop;
pub mod moleprop;
pub mod orbtprop;

#[pymodule]
fn core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    base::register_module(m)?;
    space::register_module(m)?;
    reader::register_module(m)?;
    utils::register_module(m)?;
    march::register_module(m)?;
    maths::register_module(m)?;
    integ::register_module(m)?;
    matrix::register_module(m)?;
    tools::register_module(m)?;

    atomprop::register_module(m)?;
    bondprop::register_module(m)?;
    fragprop::register_module(m)?;
    gridprop::register_module(m)?;
    moleprop::register_module(m)?;
    orbtprop::register_module(m)?;
    Ok(())
}
