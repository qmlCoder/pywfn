use numpy::{IntoPyArray, PyArray2};
use pyo3::prelude::*;
use rswfn::{self, datas::elems};

#[pyclass]
pub struct ELEM {
    core: rswfn::datas::elems::ELEM,
}

#[pymethods]
impl ELEM {
    /// 获取元素的原子序数
    #[getter]
    fn idx(&self) -> isize {
        self.core.idx
    }

    /// 获取元素的符号（如 "H", "C", "O"）
    #[getter]
    fn sym(&self) -> String {
        self.core.sym.clone()
    }

    /// 获取元素的原子半径（单位：Å）
    #[getter]
    fn rad(&self) -> f64 {
        self.core.rad
    }

    /// 获取元素在周期表中的周期数
    #[getter]
    fn row(&self) -> usize {
        self.core.row
    }

    /// 获取元素在周期表中的族数
    #[getter]
    fn col(&self) -> usize {
        self.core.col
    }
    #[getter]
    fn color(&self) -> String {
        self.core.color.clone()
    }
}

/// 获取所有元素的列表
///
/// ## 返回
/// 返回包含所有元素的 `ELEM` 对象列表
#[pyfunction]
pub fn ELEMS() -> Vec<ELEM> {
    let elems: Vec<ELEM> = rswfn::datas::elems::ELEMS
        .clone()
        .into_iter()
        .map(|e| ELEM { core: e })
        .collect();
    elems
}

/// 通过元素符号查找元素信息
///
/// ## 参数
/// - `sym`: 元素符号（如 "H", "C", "O"）
///
/// ## 返回
/// 返回对应的 `ELEM` 对象
#[pyfunction]
pub fn get_by_sym(sym: &str) -> ELEM {
    let elem = rswfn::datas::elems::get_by_sym(sym);
    ELEM { core: elem }
}

#[pyfunction]
pub fn get_by_idx(idx: usize) -> ELEM {
    let elem = rswfn::datas::elems::get_by_idx(idx);
    ELEM { core: elem }
}

#[pyfunction]
pub fn NLMS() -> [(usize, usize, isize); 59] {
    rswfn::datas::elems::NLMS
}

#[pyfunction]
pub fn ELES() -> [[usize; 59]; 118] {
    rswfn::datas::elems::ELES
}

#[pyfunction]
pub fn get_ab_nele(atm: usize) -> [usize; 2] {
    rswfn::datas::elems::get_ab_nele(atm)
}

#[pyfunction]
pub fn OXIS() -> Vec<Vec<isize>> {
    rswfn::datas::elems::OXIS
        .iter()
        .map(|row| row.to_vec())
        .collect()
}

/// 根据元素符号返回原子序数
#[pyfunction]
pub fn sym2atm(sym: &str) -> usize {
    rswfn::datas::elems::sym2atm(sym)
}

/// 根据原子序数返回元素符号
#[pyfunction]
pub fn atm2sym(atm: usize) -> String {
    rswfn::datas::elems::atm2sym(atm)
}

/// 返回指定原子的有效核电荷数（用于 pi 电子分析）
///
/// ## 参数
/// - `atm`: 原子序号
/// - `neb`: 参与共轭的电子数
#[pyfunction]
pub fn zeff(atm: usize, neb: usize) -> f64 {
    rswfn::datas::elems::zeff(atm, neb)
}

/// 根据原子序数返回原子半径
#[pyfunction]
pub fn atm2rad(idx: usize) -> f64 {
    rswfn::datas::elems::atm2rad(idx)
}

/// 根据原子序数返回周期数
#[pyfunction]
pub fn atm2row(atm: usize) -> usize {
    rswfn::datas::elems::atm2row(atm)
}

/// 返回原始的元组形式元素数据表
#[pyfunction]
pub fn ELEM_DATA() -> Vec<(isize, String, f64, usize, usize, String)> {
    rswfn::datas::elems::ELEM_DATA
        .iter()
        .map(|&(idx, sym, rad, row, col, color)| {
            (idx, sym.to_string(), rad, row, col, color.to_string())
        })
        .collect()
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "elems")?;
    m.add_class::<ELEM>()?;
    m.add_function(wrap_pyfunction!(ELEMS, &m)?)?;
    m.add_function(wrap_pyfunction!(get_by_sym, &m)?)?;
    m.add_function(wrap_pyfunction!(get_by_idx, &m)?)?;
    m.add_function(wrap_pyfunction!(NLMS, &m)?)?;
    m.add_function(wrap_pyfunction!(ELES, &m)?)?;
    m.add_function(wrap_pyfunction!(get_ab_nele, &m)?)?;
    m.add_function(wrap_pyfunction!(OXIS, &m)?)?;
    m.add_function(wrap_pyfunction!(sym2atm, &m)?)?;
    m.add_function(wrap_pyfunction!(atm2sym, &m)?)?;
    m.add_function(wrap_pyfunction!(zeff, &m)?)?;
    m.add_function(wrap_pyfunction!(atm2rad, &m)?)?;
    m.add_function(wrap_pyfunction!(atm2row, &m)?)?;
    m.add_function(wrap_pyfunction!(ELEM_DATA, &m)?)?;
    parent_module.add_submodule(&m)
}
