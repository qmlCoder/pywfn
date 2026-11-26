use numpy::PyArrayMethods;
use numpy::{IntoPyArray, PyArray2};
use pyo3::prelude::*;

use rswfn;

/// 获取 D 轨道（笛卡尔型）的变换矩阵
#[pyfunction]
fn DMAT(py: Python) -> Py<PyArray2<f64>> {
    rswfn::datas::btrans::DMAT.clone().into_pyarray(py).unbind()
}

/// 获取 F 轨道（笛卡尔型）的变换矩阵
#[pyfunction]
fn FMAT(py: Python) -> Py<PyArray2<f64>> {
    rswfn::datas::btrans::FMAT.clone().into_pyarray(py).unbind()
}

/// 获取 G 轨道（笛卡尔型）的变换矩阵
#[pyfunction]
fn GMAT(py: Python) -> Py<PyArray2<f64>> {
    rswfn::datas::btrans::GMAT.clone().into_pyarray(py).unbind()
}

/// 获取 H 轨道（笛卡尔型）的变换矩阵
#[pyfunction]
fn HMAT(py: Python) -> Py<PyArray2<f64>> {
    rswfn::datas::btrans::HMAT.clone().into_pyarray(py).unbind()
}

#[pyfunction]
fn sph2car_tmat(py: Python, ang: usize) -> Py<PyArray2<f64>> {
    let tmat = rswfn::datas::btrans::sph2car_tmat(ang);
    tmat.into_pyarray(py).unbind()
}

#[pyfunction]
fn car2sph_tmat(py: Python, ang: usize) -> Py<PyArray2<f64>> {
    let tmat = rswfn::datas::btrans::car2sph_tmat(ang);
    tmat.into_pyarray(py).unbind()
}

#[pyfunction]
fn car2sph_fmat(
    py: Python,
    smat_car: Bound<'_, PyArray2<f64>>,
    angs: Vec<usize>,
) -> Py<PyArray2<f64>> {
    let smat_car = smat_car.readonly().as_array().to_owned();
    let fmat = rswfn::datas::btrans::car2sph_fmat(&smat_car, &angs);
    fmat.into_pyarray(py).unbind()
}

// ===== SYMS 函数 =====

/// 获取混合 S 轨道的对称性标签列表
#[pyfunction]
fn MIX_S_SYMS() -> Vec<String> {
    rswfn::datas::btrans::MIX_S_SYMS
        .iter()
        .map(|&s| s.to_string())
        .collect()
}

/// 获取混合 P 轨道的对称性标签列表
#[pyfunction]
fn MIX_P_SYMS() -> Vec<String> {
    rswfn::datas::btrans::MIX_P_SYMS
        .iter()
        .map(|&s| s.to_string())
        .collect()
}

/// 获取笛卡尔 D 轨道的对称性标签列表
#[pyfunction]
fn CAR_D_SYMS() -> Vec<String> {
    rswfn::datas::btrans::CAR_D_SYMS
        .iter()
        .map(|&s| s.to_string())
        .collect()
}

/// 获取笛卡尔 F 轨道的对称性标签列表
#[pyfunction]
fn CAR_F_SYMS() -> Vec<String> {
    rswfn::datas::btrans::CAR_F_SYMS
        .iter()
        .map(|&s| s.to_string())
        .collect()
}

/// 获取笛卡尔 G 轨道的对称性标签列表
#[pyfunction]
fn CAR_G_SYMS() -> Vec<String> {
    rswfn::datas::btrans::CAR_G_SYMS
        .iter()
        .map(|&s| s.to_string())
        .collect()
}

/// 获取笛卡尔 H 轨道的对称性标签列表
#[pyfunction]
fn CAR_H_SYMS() -> Vec<String> {
    rswfn::datas::btrans::CAR_H_SYMS
        .iter()
        .map(|&s| s.to_string())
        .collect()
}

/// 获取球谐 D 轨道的对称性标签列表
#[pyfunction]
fn SPH_D_SYMS() -> Vec<String> {
    rswfn::datas::btrans::SPH_D_SYMS
        .iter()
        .map(|&s| s.to_string())
        .collect()
}

/// 获取球谐 F 轨道的对称性标签列表
#[pyfunction]
fn SPH_F_SYMS() -> Vec<String> {
    rswfn::datas::btrans::SPH_F_SYMS
        .iter()
        .map(|&s| s.to_string())
        .collect()
}

/// 获取球谐 G 轨道的对称性标签列表
#[pyfunction]
fn SPH_G_SYMS() -> Vec<String> {
    rswfn::datas::btrans::SPH_G_SYMS
        .iter()
        .map(|&s| s.to_string())
        .collect()
}

/// 获取球谐 H 轨道的对称性标签列表
#[pyfunction]
fn SPH_H_SYMS() -> Vec<String> {
    rswfn::datas::btrans::SPH_H_SYMS
        .iter()
        .map(|&s| s.to_string())
        .collect()
}

// ===== LMNS 函数 =====

/// 获取混合 S 轨道的量子数 (l, m, n) 列表
#[pyfunction]
fn MIX_S_LMNS() -> Vec<Vec<usize>> {
    rswfn::datas::btrans::MIX_S_LMNS
        .iter()
        .map(|arr| arr.to_vec())
        .collect()
}

/// 获取混合 P 轨道的量子数 (l, m, n) 列表
#[pyfunction]
fn MIX_P_LMNS() -> Vec<Vec<usize>> {
    rswfn::datas::btrans::MIX_P_LMNS
        .iter()
        .map(|arr| arr.to_vec())
        .collect()
}

/// 获取笛卡尔 D 轨道的量子数 (l, m, n) 列表
#[pyfunction]
fn CAR_D_LMNS() -> Vec<Vec<usize>> {
    rswfn::datas::btrans::CAR_D_LMNS
        .iter()
        .map(|arr| arr.to_vec())
        .collect()
}

/// 获取笛卡尔 F 轨道的量子数 (l, m, n) 列表
#[pyfunction]
fn CAR_F_LMNS() -> Vec<Vec<usize>> {
    rswfn::datas::btrans::CAR_F_LMNS
        .iter()
        .map(|arr| arr.to_vec())
        .collect()
}

/// 获取笛卡尔 G 轨道的量子数 (l, m, n) 列表
#[pyfunction]
fn CAR_G_LMNS() -> Vec<Vec<usize>> {
    rswfn::datas::btrans::CAR_G_LMNS
        .iter()
        .map(|arr| arr.to_vec())
        .collect()
}

/// 获取笛卡尔 H 轨道的量子数 (l, m, n) 列表
#[pyfunction]
fn CAR_H_LMNS() -> Vec<Vec<usize>> {
    rswfn::datas::btrans::CAR_H_LMNS
        .iter()
        .map(|arr| arr.to_vec())
        .collect()
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "btrans")?;

    m.add_function(wrap_pyfunction!(DMAT, &m)?)?;
    m.add_function(wrap_pyfunction!(FMAT, &m)?)?;
    m.add_function(wrap_pyfunction!(GMAT, &m)?)?;
    m.add_function(wrap_pyfunction!(HMAT, &m)?)?;

    m.add_function(wrap_pyfunction!(sph2car_tmat, &m)?)?;
    m.add_function(wrap_pyfunction!(car2sph_tmat, &m)?)?;
    m.add_function(wrap_pyfunction!(car2sph_fmat, &m)?)?;

    m.add_function(wrap_pyfunction!(MIX_S_SYMS, &m)?)?;
    m.add_function(wrap_pyfunction!(MIX_P_SYMS, &m)?)?;
    m.add_function(wrap_pyfunction!(CAR_D_SYMS, &m)?)?;
    m.add_function(wrap_pyfunction!(CAR_F_SYMS, &m)?)?;
    m.add_function(wrap_pyfunction!(CAR_G_SYMS, &m)?)?;
    m.add_function(wrap_pyfunction!(CAR_H_SYMS, &m)?)?;
    m.add_function(wrap_pyfunction!(SPH_D_SYMS, &m)?)?;
    m.add_function(wrap_pyfunction!(SPH_F_SYMS, &m)?)?;
    m.add_function(wrap_pyfunction!(SPH_G_SYMS, &m)?)?;
    m.add_function(wrap_pyfunction!(SPH_H_SYMS, &m)?)?;

    m.add_function(wrap_pyfunction!(MIX_S_LMNS, &m)?)?;
    m.add_function(wrap_pyfunction!(MIX_P_LMNS, &m)?)?;
    m.add_function(wrap_pyfunction!(CAR_D_LMNS, &m)?)?;
    m.add_function(wrap_pyfunction!(CAR_F_LMNS, &m)?)?;
    m.add_function(wrap_pyfunction!(CAR_G_LMNS, &m)?)?;
    m.add_function(wrap_pyfunction!(CAR_H_LMNS, &m)?)?;

    parent_module.add_submodule(&m)
}
