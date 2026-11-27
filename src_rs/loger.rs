#![allow(non_snake_case)]

use pyo3::prelude::*;
use rswfn;

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "loger")?;
    m.add_function(wrap_pyfunction!(init_log, &m)?)?;
    m.add_function(wrap_pyfunction!(init_stderr_log, &m)?)?;
    parent_module.add_submodule(&m)
}

/// 初始化日志，输出到日志文件。
#[pyfunction]
pub fn init_log() {
    rswfn::loger::init_log()
}

/// 初始化日志，输出到标准错误流。
#[pyfunction]
pub fn init_stderr_log() {
    rswfn::loger::init_stderr_log()
}
