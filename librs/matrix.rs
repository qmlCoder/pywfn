use pyo3::prelude::*;
use rswfn;

#[pyfunction]
pub fn ele_mat(matc: Vec<Vec<f64>>, mats: Vec<Vec<f64>>) -> PyResult<Vec<Vec<f64>>> {
    Ok(rswfn::matrix::ele_mat(&matc, &mats))
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "matrix")?; // 为当前的rs文件创建一个子模块
    m.add_function(wrap_pyfunction!(ele_mat, &m)?)?; // 添加这个rs文件中的函数到子模块中
    parent_module.add_submodule(&m) // 将子模块添加到父模块中
}
