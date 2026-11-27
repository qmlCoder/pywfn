use numpy::{IntoPyArray, PyArray2, PyArrayMethods};
use pyo3::prelude::*;
use rswfn;

use crate::base::Mole;

/// 分子轨道系数的用途
#[pyclass(from_py_object)]
#[derive(Clone)]
pub struct MocvDest {
    pub kind: MocvKind,
}

/// 分子轨道系数用途的内部表示
#[derive(Clone)]
pub enum MocvKind {
    PiProperty,
    OrderDecom(OrderDType),
    RoundOrder,
}

impl MocvDest {
    /// 转换为 rswfn 中的枚举
    pub fn to_core(&self) -> rswfn::props::orbtprop::mocv::MocvDest {
        match &self.kind {
            MocvKind::PiProperty => rswfn::props::orbtprop::mocv::MocvDest::PiProperty,
            MocvKind::OrderDecom(dtype) => {
                rswfn::props::orbtprop::mocv::MocvDest::OrderDecom(dtype.to_core())
            }
            MocvKind::RoundOrder => rswfn::props::orbtprop::mocv::MocvDest::RoundOrder,
        }
    }
}

#[pymethods]
impl MocvDest {
    /// 用于计算 π 电子相关的性质
    #[staticmethod]
    pub fn PiProperty() -> Self {
        Self {
            kind: MocvKind::PiProperty,
        }
    }

    /// 用于分解键级
    ///
    /// ## 参数
    /// - `dtype`: 需要分解的键级类型
    #[staticmethod]
    pub fn OrderDecom(dtype: OrderDType) -> Self {
        Self {
            kind: MocvKind::OrderDecom(dtype),
        }
    }

    /// 用于取整键级
    #[staticmethod]
    pub fn RoundOrder() -> Self {
        Self {
            kind: MocvKind::RoundOrder,
        }
    }
}

/// 键级分解的类型
#[pyclass(eq, eq_int, from_py_object)]
#[derive(Clone, PartialEq)]
pub enum OrderDType {
    /// 总键级
    SI,
    /// π 键级
    PI,
    /// σ 键级
    DE,
    /// 相位键级
    PH,
}

impl OrderDType {
    /// 转换为 rswfn 中的枚举
    pub fn to_core(&self) -> rswfn::props::orbtprop::mocv::OrderDType {
        match self {
            OrderDType::SI => rswfn::props::orbtprop::mocv::OrderDType::SI,
            OrderDType::PI => rswfn::props::orbtprop::mocv::OrderDType::PI,
            OrderDType::DE => rswfn::props::orbtprop::mocv::OrderDType::DE,
            OrderDType::PH => rswfn::props::orbtprop::mocv::OrderDType::PH,
        }
    }
}

/// 分子轨道系数变换矩阵
///
/// ## 参数
/// - `py`: Python 解释器标记
/// - `L`: 分子轨道系数矩阵
/// - `rcoefs`: 参考系数矩阵
/// - `keeps`: 保留的轨道掩码
#[pyfunction]
fn mocv(
    py: Python,
    L: Bound<'_, PyArray2<f64>>,
    rcoefs: Bound<'_, PyArray2<f64>>,
    keeps: Bound<'_, PyArray2<f64>>,
) -> Py<PyArray2<f64>> {
    let L = L.readonly().as_array().to_owned();
    let rcoefs = rcoefs.readonly().as_array().to_owned();
    let keeps = keeps.readonly().as_array().to_owned();
    rswfn::props::orbtprop::mocv::mocv(&L, &rcoefs, &keeps)
        .into_pyarray(py)
        .unbind()
}

/// 获取参与分子轨道系数变换的原子索引
///
/// ## 参数
/// - `mole`: 分子对象
/// - `dest`: 分子轨道系数的用途
#[pyfunction]
fn mocv_atos(mole: Mole, dest: MocvDest) -> Vec<usize> {
    rswfn::props::orbtprop::mocv::mocv_atos(&mole.core, &dest.to_core())
}

/// 获取参与分子轨道系数变换的基函数符号
///
/// ## 参数
/// - `dest`: 分子轨道系数的用途
/// - `geom`: 几何结构类型
#[pyfunction]
fn mocv_syms(dest: MocvDest, geom: usize) -> Vec<String> {
    rswfn::props::orbtprop::mocv::mocv_syms(&dest.to_core(), geom)
}

#[pyfunction]
fn get_tmat(py: Python, num: usize, L: Bound<'_, PyArray2<f64>>) -> Py<PyArray2<f64>> {
    let L = L.readonly().as_array().to_owned();
    let tmat = rswfn::props::orbtprop::mocv::get_tmat(num, &L);
    tmat.into_pyarray(py).unbind()
}

#[pyfunction]
fn get_tmat_sph(py: Python, num: usize, L: Bound<'_, PyArray2<f64>>) -> Py<PyArray2<f64>> {
    let L = L.readonly().as_array().to_owned();
    let tmat = rswfn::props::orbtprop::mocv::get_tmat_sph(num, &L);
    tmat.into_pyarray(py).unbind()
}

#[pyfunction]
fn get_tmat_car(py: Python, num: usize, L: Bound<'_, PyArray2<f64>>) -> Py<PyArray2<f64>> {
    let L = L.readonly().as_array().to_owned();
    let tmat = rswfn::props::orbtprop::mocv::get_tmat_car(num, &L);
    tmat.into_pyarray(py).unbind()
}

#[pyfunction]
fn tmat_s(py: Python, L: Bound<'_, PyArray2<f64>>) -> Py<PyArray2<f64>> {
    let L = L.readonly().as_array().to_owned();
    let tmat = rswfn::props::orbtprop::mocv::tmat_s(&L);
    tmat.into_pyarray(py).unbind()
}

#[pyfunction]
fn tmat_p(py: Python, L: Bound<'_, PyArray2<f64>>) -> Py<PyArray2<f64>> {
    let L = L.readonly().as_array().to_owned();
    let tmat = rswfn::props::orbtprop::mocv::tmat_p(&L);
    tmat.into_pyarray(py).unbind()
}

#[pyfunction]
fn tmat_d(py: Python, L: Bound<'_, PyArray2<f64>>) -> Py<PyArray2<f64>> {
    let L = L.readonly().as_array().to_owned();
    let tmat = rswfn::props::orbtprop::mocv::tmat_d(&L);
    tmat.into_pyarray(py).unbind()
}

#[pyfunction]
fn tmat_f(py: Python, L: Bound<'_, PyArray2<f64>>) -> Py<PyArray2<f64>> {
    let L = L.readonly().as_array().to_owned();
    let tmat = rswfn::props::orbtprop::mocv::tmat_f(&L);
    tmat.into_pyarray(py).unbind()
}

#[pyfunction]
fn tmat_g(py: Python, L: Bound<'_, PyArray2<f64>>) -> Py<PyArray2<f64>> {
    let L = L.readonly().as_array().to_owned();
    let tmat = rswfn::props::orbtprop::mocv::tmat_g(&L);
    tmat.into_pyarray(py).unbind()
}

#[pyfunction]
fn get_nmat(py: Python, lmns: Vec<[usize; 3]>) -> Py<PyArray2<f64>> {
    let nmat = rswfn::props::orbtprop::mocv::get_nmat(&lmns);
    nmat.into_pyarray(py).unbind()
}

#[pyfunction]
fn nmat_d(py: Python) -> Py<PyArray2<f64>> {
    let nmat = rswfn::props::orbtprop::mocv::nmat_d();
    nmat.into_pyarray(py).unbind()
}

#[pyfunction]
fn nmat_f(py: Python) -> Py<PyArray2<f64>> {
    let nmat = rswfn::props::orbtprop::mocv::nmat_f();
    nmat.into_pyarray(py).unbind()
}

#[pyfunction]
fn nmat_g(py: Python) -> Py<PyArray2<f64>> {
    let nmat = rswfn::props::orbtprop::mocv::nmat_g();
    nmat.into_pyarray(py).unbind()
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "mocv")?;
    m.add_class::<MocvDest>()?;
    m.add_class::<OrderDType>()?;
    m.add_function(wrap_pyfunction!(mocv, &m)?)?;
    m.add_function(wrap_pyfunction!(mocv_atos, &m)?)?;
    m.add_function(wrap_pyfunction!(mocv_syms, &m)?)?;
    m.add_function(wrap_pyfunction!(get_tmat, &m)?)?;
    m.add_function(wrap_pyfunction!(get_tmat_sph, &m)?)?;
    m.add_function(wrap_pyfunction!(get_tmat_car, &m)?)?;
    m.add_function(wrap_pyfunction!(tmat_s, &m)?)?;
    m.add_function(wrap_pyfunction!(tmat_p, &m)?)?;
    m.add_function(wrap_pyfunction!(tmat_d, &m)?)?;
    m.add_function(wrap_pyfunction!(tmat_f, &m)?)?;
    m.add_function(wrap_pyfunction!(tmat_g, &m)?)?;
    m.add_function(wrap_pyfunction!(get_nmat, &m)?)?;
    m.add_function(wrap_pyfunction!(nmat_d, &m)?)?;
    m.add_function(wrap_pyfunction!(nmat_f, &m)?)?;
    m.add_function(wrap_pyfunction!(nmat_g, &m)?)?;
    parent_module.add_submodule(&m)
}
