use numpy::{IntoPyArray, PyArray2, PyArrayMethods};
use pyo3::prelude::*;
use rswfn;

#[pyfunction]
fn get_tmat(py: Python, num: usize, L: Bound<'_, PyArray2<f64>>) -> Py<PyArray2<f64>> {
    let L = L.readonly().as_array().to_owned();
    let tmat = rswfn::datas::mocv::get_tmat(num, &L);
    tmat.into_pyarray(py).unbind()
}

#[pyfunction]
fn get_tmat_sph(py: Python, num: usize, L: Bound<'_, PyArray2<f64>>) -> Py<PyArray2<f64>> {
    let L = L.readonly().as_array().to_owned();
    let tmat = rswfn::datas::mocv::get_tmat_sph(num, &L);
    tmat.into_pyarray(py).unbind()
}

#[pyfunction]
fn get_tmat_car(py: Python, num: usize, L: Bound<'_, PyArray2<f64>>) -> Py<PyArray2<f64>> {
    let L = L.readonly().as_array().to_owned();
    let tmat = rswfn::datas::mocv::get_tmat_car(num, &L);
    tmat.into_pyarray(py).unbind()
}

#[pyfunction]
fn tmat_s(py: Python, L: Bound<'_, PyArray2<f64>>) -> Py<PyArray2<f64>> {
    let L = L.readonly().as_array().to_owned();
    let tmat = rswfn::datas::mocv::tmat_s(&L);
    tmat.into_pyarray(py).unbind()
}

#[pyfunction]
fn tmat_p(py: Python, L: Bound<'_, PyArray2<f64>>) -> Py<PyArray2<f64>> {
    let L = L.readonly().as_array().to_owned();
    let tmat = rswfn::datas::mocv::tmat_p(&L);
    tmat.into_pyarray(py).unbind()
}

#[pyfunction]
fn tmat_d(py: Python, L: Bound<'_, PyArray2<f64>>) -> Py<PyArray2<f64>> {
    let L = L.readonly().as_array().to_owned();
    let tmat = rswfn::datas::mocv::tmat_d(&L);
    tmat.into_pyarray(py).unbind()
}

#[pyfunction]
fn tmat_f(py: Python, L: Bound<'_, PyArray2<f64>>) -> Py<PyArray2<f64>> {
    let L = L.readonly().as_array().to_owned();
    let tmat = rswfn::datas::mocv::tmat_f(&L);
    tmat.into_pyarray(py).unbind()
}

#[pyfunction]
fn tmat_g(py: Python, L: Bound<'_, PyArray2<f64>>) -> Py<PyArray2<f64>> {
    let L = L.readonly().as_array().to_owned();
    let tmat = rswfn::datas::mocv::tmat_g(&L);
    tmat.into_pyarray(py).unbind()
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "mocv")?;
    m.add_function(wrap_pyfunction!(get_tmat, &m)?)?;
    m.add_function(wrap_pyfunction!(get_tmat_sph, &m)?)?;
    m.add_function(wrap_pyfunction!(get_tmat_car, &m)?)?;
    m.add_function(wrap_pyfunction!(tmat_s, &m)?)?;
    m.add_function(wrap_pyfunction!(tmat_p, &m)?)?;
    m.add_function(wrap_pyfunction!(tmat_d, &m)?)?;
    m.add_function(wrap_pyfunction!(tmat_f, &m)?)?;
    m.add_function(wrap_pyfunction!(tmat_g, &m)?)?;
    parent_module.add_submodule(&m)
}
