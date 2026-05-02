use numpy::{IntoPyArray, PyArray2};
use pyo3::prelude::*;

use rswfn;

#[pyfunction]
fn DMAT(py: Python) -> Py<PyArray2<f64>> {
    rswfn::datas::btrans::DMAT.clone().into_pyarray(py).unbind()
}

#[pyfunction]
fn FMAT(py: Python) -> Py<PyArray2<f64>> {
    rswfn::datas::btrans::FMAT.clone().into_pyarray(py).unbind()
}

#[pyfunction]
fn GMAT(py: Python) -> Py<PyArray2<f64>> {
    rswfn::datas::btrans::GMAT.clone().into_pyarray(py).unbind()
}

#[pyfunction]
fn HMAT(py: Python) -> Py<PyArray2<f64>> {
    rswfn::datas::btrans::HMAT.clone().into_pyarray(py).unbind()
}

#[pyfunction]
fn MIX_S_LMNS() -> PyResult<Vec<Vec<usize>>> {
    // 将数组转换为 Vec<Vec<usize>> 返回给 Python
    let result: Vec<Vec<usize>> = rswfn::datas::btrans::MIX_S_LMNS
        .iter()
        .map(|arr| arr.to_vec())
        .collect();
    Ok(result)
}

#[pyfunction]
fn MIX_P_LMNS() -> PyResult<Vec<Vec<usize>>> {
    // 将数组转换为 Vec<Vec<usize>> 返回给 Python
    let result: Vec<Vec<usize>> = rswfn::datas::btrans::MIX_P_LMNS
        .iter()
        .map(|arr| arr.to_vec())
        .collect();
    Ok(result)
}

#[pyfunction]
fn CAR_D_LMNS() -> PyResult<Vec<Vec<usize>>> {
    // 将数组转换为 Vec<Vec<usize>> 返回给 Python
    let result: Vec<Vec<usize>> = rswfn::datas::btrans::CAR_D_LMNS
        .iter()
        .map(|arr| arr.to_vec())
        .collect();
    Ok(result)
}

#[pyfunction]
fn CAR_F_LMNS() -> PyResult<Vec<Vec<usize>>> {
    // 将数组转换为 Vec<Vec<usize>> 返回给 Python
    let result: Vec<Vec<usize>> = rswfn::datas::btrans::CAR_F_LMNS
        .iter()
        .map(|arr| arr.to_vec())
        .collect();
    Ok(result)
}

#[pyfunction]
fn CAR_G_LMNS() -> PyResult<Vec<Vec<usize>>> {
    // 将数组转换为 Vec<Vec<usize>> 返回给 Python
    let result: Vec<Vec<usize>> = rswfn::datas::btrans::CAR_G_LMNS
        .iter()
        .map(|arr| arr.to_vec())
        .collect();
    Ok(result)
}

#[pyfunction]
fn CAR_H_LMNS() -> PyResult<Vec<Vec<usize>>> {
    // 将数组转换为 Vec<Vec<usize>> 返回给 Python
    let result: Vec<Vec<usize>> = rswfn::datas::btrans::CAR_H_LMNS
        .iter()
        .map(|arr| arr.to_vec())
        .collect();
    Ok(result)
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "btrans")?;
    m.add_function(wrap_pyfunction!(DMAT, &m)?)?;
    m.add_function(wrap_pyfunction!(FMAT, &m)?)?;
    m.add_function(wrap_pyfunction!(GMAT, &m)?)?;
    m.add_function(wrap_pyfunction!(HMAT, &m)?)?;

    m.add_function(wrap_pyfunction!(MIX_S_LMNS, &m)?)?;
    m.add_function(wrap_pyfunction!(MIX_P_LMNS, &m)?)?;
    m.add_function(wrap_pyfunction!(CAR_D_LMNS, &m)?)?;
    m.add_function(wrap_pyfunction!(CAR_F_LMNS, &m)?)?;
    m.add_function(wrap_pyfunction!(CAR_G_LMNS, &m)?)?;
    m.add_function(wrap_pyfunction!(CAR_H_LMNS, &m)?)?;
    parent_module.add_submodule(&m)
}
