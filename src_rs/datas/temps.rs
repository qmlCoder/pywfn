use pyo3::prelude::*;
use rswfn;

/// Gaussian 输入文件模板
#[pyfunction]
fn GJF() -> String {
    rswfn::datas::temps::GJF.to_string()
}

/// 补充材料（SI）文档模板
#[pyfunction]
fn SI() -> String {
    rswfn::datas::temps::SI.to_string()
}

/// XYZ 文件模板
#[pyfunction]
fn XYZ() -> String {
    rswfn::datas::temps::XYZ.to_string()
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "temps")?;
    m.add_function(wrap_pyfunction!(GJF, &m)?)?;
    m.add_function(wrap_pyfunction!(SI, &m)?)?;
    m.add_function(wrap_pyfunction!(XYZ, &m)?)?;
    parent_module.add_submodule(&m)
}
