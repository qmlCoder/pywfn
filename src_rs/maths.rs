use numpy::{IntoPyArray, PyArray2, PyArrayMethods};
use pyo3::prelude::*;
use rswfn;

#[pyfunction]
pub fn lag_intpol(xs: Vec<f64>, ys: Vec<f64>, ts: Vec<f64>) -> Vec<f64> {
    ts.iter()
        .map(|t| rswfn::maths::lag_intpol(&xs, &ys, *t))
        .collect()
}

#[pyfunction]
pub fn matrix_phalf(py: Python, mat: Bound<'_, PyArray2<f64>>) -> Py<PyArray2<f64>> {
    let mat = mat.readonly().as_array().to_owned();
    let phalf = rswfn::maths::matrix_phalf(&mat);
    phalf.into_pyarray(py).unbind()
}

#[pyfunction]
pub fn matrix_nhalf(py: Python, mat: Bound<'_, PyArray2<f64>>) -> Py<PyArray2<f64>> {
    let mat = mat.readonly().as_array().to_owned();
    let nhalf = rswfn::maths::matrix_nhalf(&mat);
    nhalf.into_pyarray(py).unbind()
}

#[pyfunction]
pub fn calc_dmat(
    py: Python,
    cmat: Bound<'_, PyArray2<f64>>,
    obt_occs: Vec<f64>,
) -> Py<PyArray2<f64>> {
    let cmat = cmat.readonly().as_array().to_owned();
    let dmat = rswfn::maths::calc_dmat(&cmat, &obt_occs);
    dmat.into_pyarray(py).unbind()
}

#[pyfunction]
pub fn fac(n: isize) -> isize {
    rswfn::maths::fac(n)
}

#[pyfunction]
pub fn fac2(n: isize) -> isize {
    rswfn::maths::fac2(n)
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "maths")?; // 为当前的rs文件创建一个子模块
    m.add_function(wrap_pyfunction!(lag_intpol, &m)?)?; //
    m.add_function(wrap_pyfunction!(calc_dmat, &m)?)?; //
    m.add_function(wrap_pyfunction!(matrix_phalf, &m)?)?; //
    m.add_function(wrap_pyfunction!(matrix_nhalf, &m)?)?; //
    m.add_function(wrap_pyfunction!(fac, &m)?)?; //
    m.add_function(wrap_pyfunction!(fac2, &m)?)?; //
    parent_module.add_submodule(&m) // 将子模块添加到父模块中
}
