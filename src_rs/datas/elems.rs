use numpy::{IntoPyArray, PyArray2};
use pyo3::prelude::*;
use rswfn::{self, datas::elems};

#[pyclass]
pub struct ELEM {
    core: rswfn::datas::elems::ELEM,
}

#[pymethods]
impl ELEM {
    #[getter]
    fn idx(&self) -> usize {
        self.core.idx
    }

    #[getter]
    fn sym(&self) -> String {
        self.core.sym.clone()
    }

    #[getter]
    fn rad(&self) -> f64 {
        self.core.rad
    }

    #[getter]
    fn row(&self) -> usize {
        self.core.row
    }

    #[getter]
    fn col(&self) -> usize {
        self.core.col
    }
}

#[pyfunction]
pub fn ELEMS() -> Vec<ELEM> {
    let elems: Vec<ELEM> = rswfn::datas::elems::ELEMS
        .clone()
        .into_iter()
        .map(|e| ELEM { core: e })
        .collect();
    elems
}

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

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "elems")?;
    m.add_class::<ELEM>()?;
    m.add_function(wrap_pyfunction!(ELEMS, &m)?)?;
    m.add_function(wrap_pyfunction!(get_by_sym, &m)?)?;
    m.add_function(wrap_pyfunction!(get_by_idx, &m)?)?;
    parent_module.add_submodule(&m)
}
