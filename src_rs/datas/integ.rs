use numpy::{IntoPyArray, PyArray2};
use pyo3::prelude::*;
use rswfn;

/// 获取径向哈密顿矩阵
#[pyfunction]
fn RHM(py: Python) -> Py<PyArray2<f64>> {
    rswfn::datas::integ::RHM.clone().into_pyarray(py).unbind()
}

/// 获取径向权重矩阵
#[pyfunction]
fn WHM(py: Python) -> Py<PyArray2<f64>> {
    rswfn::datas::integ::WHM.clone().into_pyarray(py).unbind()
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "integ")?;
    m.add_function(wrap_pyfunction!(RHM, &m)?)?;
    m.add_function(wrap_pyfunction!(WHM, &m)?)?;
    parent_module.add_submodule(&m)
}
