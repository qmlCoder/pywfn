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

// ===== 全局符号表 =====

/// 获取所有基函数的对称性标签列表
#[pyfunction]
fn SYMS() -> Vec<String> {
    rswfn::datas::btrans::SYMS
        .iter()
        .map(|&s| s.to_string())
        .collect()
}

/// 获取所有基函数的量子数 (l, m, n) 列表
#[pyfunction]
fn LMNS() -> Vec<Vec<usize>> {
    rswfn::datas::btrans::LMNS
        .iter()
        .map(|arr| arr.to_vec())
        .collect()
}

/// 获取指定角动量的量子数 (l, m, n) 列表
///
/// ## 参数
/// - `ang`: 角动量量子数
#[pyfunction]
fn ang2lmns(ang: usize) -> Vec<Vec<usize>> {
    rswfn::datas::btrans::ang2lmns(ang)
        .into_iter()
        .map(|arr| arr.to_vec())
        .collect()
}

/// 获取指定角动量的对称性标签列表
///
/// ## 参数
/// - `ang`: 角动量量子数
/// - `form`: 基函数形式（"sph" 或 "car"）
#[pyfunction]
fn ang2syms(ang: usize, form: &str) -> Vec<String> {
    rswfn::datas::btrans::ang2syms(ang, form)
        .into_iter()
        .map(|s| s.to_string())
        .collect()
}

/// 获取笛卡尔型基函数的对称性标签（带缓存）
#[pyfunction]
fn get_car_syms() -> Vec<String> {
    rswfn::datas::btrans::get_car_syms()
        .into_iter()
        .map(|s| s.to_string())
        .collect()
}

/// 获取笛卡尔型基函数的对称性标签（不使用缓存）
#[pyfunction]
fn get_car_syms_no_cache() -> Vec<String> {
    rswfn::datas::btrans::get_car_syms_no_cache()
        .into_iter()
        .map(|s| s.to_string())
        .collect()
}

/// 获取笛卡尔型基函数的对称性标签（预填充缓存）
#[pyfunction]
fn get_car_syms_prime_cache() -> Vec<String> {
    rswfn::datas::btrans::get_car_syms_prime_cache()
        .into_iter()
        .map(|s| s.to_string())
        .collect()
}

/// 获取球谐型基函数的对称性标签（带缓存）
#[pyfunction]
fn get_sph_syms() -> Vec<String> {
    rswfn::datas::btrans::get_sph_syms()
        .into_iter()
        .map(|s| s.to_string())
        .collect()
}

/// 获取球谐型基函数的对称性标签（不使用缓存）
#[pyfunction]
fn get_sph_syms_no_cache() -> Vec<String> {
    rswfn::datas::btrans::get_sph_syms_no_cache()
        .into_iter()
        .map(|s| s.to_string())
        .collect()
}

/// 获取球谐型基函数的对称性标签（预填充缓存）
#[pyfunction]
fn get_sph_syms_prime_cache() -> Vec<String> {
    rswfn::datas::btrans::get_sph_syms_prime_cache()
        .into_iter()
        .map(|s| s.to_string())
        .collect()
}

/// 获取混合型基函数的对称性标签（带缓存）
#[pyfunction]
fn get_mix_syms() -> Vec<String> {
    rswfn::datas::btrans::get_mix_syms()
        .into_iter()
        .map(|s| s.to_string())
        .collect()
}

/// 获取混合型基函数的对称性标签（不使用缓存）
#[pyfunction]
fn get_mix_syms_no_cache() -> Vec<String> {
    rswfn::datas::btrans::get_mix_syms_no_cache()
        .into_iter()
        .map(|s| s.to_string())
        .collect()
}

/// 获取混合型基函数的对称性标签（预填充缓存）
#[pyfunction]
fn get_mix_syms_prime_cache() -> Vec<String> {
    rswfn::datas::btrans::get_mix_syms_prime_cache()
        .into_iter()
        .map(|s| s.to_string())
        .collect()
}

/// 根据基函数个数获取对应的对称性标签列表
///
/// ## 参数
/// - `len`: 基函数个数
#[pyfunction]
fn len2syms(len: usize) -> Vec<String> {
    rswfn::datas::btrans::len2syms(len)
}

/// 根据基函数个数获取对称性标签（不使用缓存）
#[pyfunction]
fn len2syms_no_cache(len: usize) -> Vec<String> {
    rswfn::datas::btrans::len2syms_no_cache(len)
}

/// 根据基函数个数获取对称性标签（预填充缓存）
#[pyfunction]
fn len2syms_prime_cache(len: usize) -> Vec<String> {
    rswfn::datas::btrans::len2syms_prime_cache(len)
}

/// 由量子数 (l, m, n) 获取对称性标签
///
/// ## 参数
/// - `lmn`: 量子数三元组 [l, m, n]
#[pyfunction]
fn lmn2sym(lmn: [usize; 3]) -> String {
    rswfn::datas::btrans::lmn2sym(&lmn)
}

/// 由对称性标签获取角动量量子数
///
/// ## 参数
/// - `sym`: 对称性标签，如 "px"
#[pyfunction]
fn sym2ang(sym: &str) -> usize {
    rswfn::datas::btrans::sym2ang(sym)
}

/// 由对称性标签获取角动量量子数（不使用缓存）
#[pyfunction]
fn sym2ang_no_cache(sym: &str) -> usize {
    rswfn::datas::btrans::sym2ang_no_cache(sym)
}

/// 由对称性标签获取角动量量子数（预填充缓存）
#[pyfunction]
fn sym2ang_prime_cache(sym: &str) -> usize {
    rswfn::datas::btrans::sym2ang_prime_cache(sym)
}

/// 由对称性标签获取量子数 (l, m, n)
///
/// ## 参数
/// - `sym`: 对称性标签，如 "px"
#[pyfunction]
fn sym2lmn(sym: &str) -> [usize; 3] {
    rswfn::datas::btrans::sym2lmn(sym)
}

/// 由对称性标签获取量子数（不使用缓存）
#[pyfunction]
fn sym2lmn_no_cache(sym: &str) -> [usize; 3] {
    rswfn::datas::btrans::sym2lmn_no_cache(sym)
}

/// 由对称性标签获取量子数（预填充缓存）
#[pyfunction]
fn sym2lmn_prime_cache(sym: &str) -> [usize; 3] {
    rswfn::datas::btrans::sym2lmn_prime_cache(sym)
}

/// 由对称性标签获取起始索引
///
/// ## 参数
/// - `sym`: 对称性标签
#[pyfunction]
fn sym2stp(sym: &str) -> usize {
    rswfn::datas::btrans::sym2stp(sym)
}

/// 按给定的基函数形式做基函数变换
///
/// ## 参数
/// - `from_smat`: 原始矩阵
/// - `from_forms`: 原始矩阵每行对应的基函数形式
/// - `into_forms`: 目标矩阵每行对应的基函数形式
/// - `angs`: 每行对应的角动量量子数
#[pyfunction]
fn get_fmat(
    py: Python,
    from_smat: Bound<'_, PyArray2<f64>>,
    from_forms: Vec<String>,
    into_forms: Vec<String>,
    angs: Vec<usize>,
) -> Py<PyArray2<f64>> {
    let from_smat = from_smat.readonly().as_array().to_owned();
    let from_forms = from_forms.iter().map(|s| s.as_str()).collect::<Vec<_>>();
    let into_forms = into_forms.iter().map(|s| s.as_str()).collect::<Vec<_>>();
    let fmat = rswfn::datas::btrans::get_fmat(&from_smat, &from_forms, &into_forms, &angs);
    fmat.into_pyarray(py).unbind()
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

    m.add_function(wrap_pyfunction!(SYMS, &m)?)?;
    m.add_function(wrap_pyfunction!(LMNS, &m)?)?;
    m.add_function(wrap_pyfunction!(ang2lmns, &m)?)?;
    m.add_function(wrap_pyfunction!(ang2syms, &m)?)?;
    m.add_function(wrap_pyfunction!(get_car_syms, &m)?)?;
    m.add_function(wrap_pyfunction!(get_car_syms_no_cache, &m)?)?;
    m.add_function(wrap_pyfunction!(get_car_syms_prime_cache, &m)?)?;
    m.add_function(wrap_pyfunction!(get_sph_syms, &m)?)?;
    m.add_function(wrap_pyfunction!(get_sph_syms_no_cache, &m)?)?;
    m.add_function(wrap_pyfunction!(get_sph_syms_prime_cache, &m)?)?;
    m.add_function(wrap_pyfunction!(get_mix_syms, &m)?)?;
    m.add_function(wrap_pyfunction!(get_mix_syms_no_cache, &m)?)?;
    m.add_function(wrap_pyfunction!(get_mix_syms_prime_cache, &m)?)?;
    m.add_function(wrap_pyfunction!(len2syms, &m)?)?;
    m.add_function(wrap_pyfunction!(len2syms_no_cache, &m)?)?;
    m.add_function(wrap_pyfunction!(len2syms_prime_cache, &m)?)?;
    m.add_function(wrap_pyfunction!(lmn2sym, &m)?)?;
    m.add_function(wrap_pyfunction!(sym2ang, &m)?)?;
    m.add_function(wrap_pyfunction!(sym2ang_no_cache, &m)?)?;
    m.add_function(wrap_pyfunction!(sym2ang_prime_cache, &m)?)?;
    m.add_function(wrap_pyfunction!(sym2lmn, &m)?)?;
    m.add_function(wrap_pyfunction!(sym2lmn_no_cache, &m)?)?;
    m.add_function(wrap_pyfunction!(sym2lmn_prime_cache, &m)?)?;
    m.add_function(wrap_pyfunction!(sym2stp, &m)?)?;
    m.add_function(wrap_pyfunction!(get_fmat, &m)?)?;

    parent_module.add_submodule(&m)
}
