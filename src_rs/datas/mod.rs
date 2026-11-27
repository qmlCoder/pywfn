mod btrans;
mod consts;
mod dhmo2;
mod elems;
mod integ;
mod lebedev;
mod rdens;
mod roots;
mod rys_root;
mod temps;

use pyo3::prelude::*;
use rswfn;

/// 阶乘相关的常数表
#[pyfunction]
fn FACS() -> Vec<f64> {
    rswfn::datas::FACS.to_vec()
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "datas")?;
    m.add_function(wrap_pyfunction!(FACS, &m)?)?;
    btrans::register_module(&m)?;
    consts::register_module(&m)?;
    dhmo2::register_module(&m)?;
    elems::register_module(&m)?;
    integ::register_module(&m)?;
    lebedev::register_module(&m)?;
    rdens::register_module(&m)?;
    roots::register_module(&m)?;
    rys_root::register_module(&m)?;
    temps::register_module(&m)?;
    parent_module.add_submodule(&m)
}
