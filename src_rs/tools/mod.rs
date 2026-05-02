use pyo3::prelude::*;

pub mod log;

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "tools")?; // 为当前的rs文件创建一个子模块
    log::register_module(&m)?;
    parent_module.add_submodule(&m) // 将子模块添加到父模块中
}
