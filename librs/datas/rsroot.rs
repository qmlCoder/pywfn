use pyo3::prelude::*;
use rswfn;

#[pyfunction]
pub fn rys_root1(X: f64) -> PyResult<(Vec<f64>, Vec<f64>)> {
    Ok(rswfn::datas::rys_root::rys_root1(X))
}

#[pyfunction]
pub fn rys_root2(X: f64) -> PyResult<(Vec<f64>, Vec<f64>)> {
    Ok(rswfn::datas::rys_root::rys_root2(X))
}

#[pyfunction]
pub fn rys_root3(X: f64) -> PyResult<(Vec<f64>, Vec<f64>)> {
    Ok(rswfn::datas::rys_root::rys_root3(X))
}

#[pyfunction]
pub fn rys_root4(X: f64) -> PyResult<(Vec<f64>, Vec<f64>)> {
    Ok(rswfn::datas::rys_root::rys_root4(X))
}

#[pyfunction]
pub fn rys_root5(X: f64) -> PyResult<(Vec<f64>, Vec<f64>)> {
    Ok(rswfn::datas::rys_root::rys_root5(X))
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "rsroot")?;
    m.add_function(wrap_pyfunction!(rys_root1, &m)?)?;
    m.add_function(wrap_pyfunction!(rys_root2, &m)?)?;
    m.add_function(wrap_pyfunction!(rys_root3, &m)?)?;
    m.add_function(wrap_pyfunction!(rys_root4, &m)?)?;
    m.add_function(wrap_pyfunction!(rys_root5, &m)?)?;
    parent_module.add_submodule(&m)
}
