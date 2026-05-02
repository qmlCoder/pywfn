use pyo3::prelude::*;
use rswfn;

#[pyfunction]
pub fn lag_intpol(xs: Vec<f64>, ys: Vec<f64>, ts: Vec<f64>) -> Vec<f64> {
    ts.iter()
        .map(|t| rswfn::maths::lag_intpol(&xs, &ys, *t))
        .collect()
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "maths")?; // 为当前的rs文件创建一个子模块
    m.add_function(wrap_pyfunction!(lag_intpol, &m)?)?; // 添加这个rs文件中的函数到子模块中
    parent_module.add_submodule(&m) // 将子模块添加到父模块中
}
