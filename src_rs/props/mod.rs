//! 各类性质（电荷、键级、能量、格点、轨道等）计算模块。

use pyo3::prelude::*;

pub mod atomprop;
pub mod bondprop;
pub mod fragprop;
pub mod gridprop;
pub mod moleprop;
pub mod orbtprop;

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "props")?;
    atomprop::register_module(&m)?;
    bondprop::register_module(&m)?;
    fragprop::register_module(&m)?;
    gridprop::register_module(&m)?;
    moleprop::register_module(&m)?;
    orbtprop::register_module(&m)?;
    parent_module.add_submodule(&m)
}
