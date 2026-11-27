#![allow(non_snake_case)]

use pyo3::prelude::*;

pub mod base;
pub mod datas;
pub mod editor;
pub mod loger;
pub mod maths;
pub mod props;
pub mod qmcal;
pub mod reader;
pub mod rsbse;
pub mod rsint;
pub mod rsxcs;
pub mod space;
pub mod tools;
pub mod utils;
pub mod writer;

#[pymodule]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    base::register_module(m)?;
    datas::register_module(m)?;
    editor::register_module(m)?;
    loger::register_module(m)?;
    maths::register_module(m)?;
    props::register_module(m)?;
    qmcal::register_module(m)?;
    reader::register_module(m)?;
    rsbse::register_module(m)?;
    rsint::register_module(m)?;
    rsxcs::register_module(m)?;
    space::register_module(m)?;
    tools::register_module(m)?;
    utils::register_module(m)?;
    writer::register_module(m)?;
    Ok(())
}
