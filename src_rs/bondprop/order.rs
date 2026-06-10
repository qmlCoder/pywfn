use pyo3::prelude::*;
use std::collections::HashMap;

use rswfn;
use rswfn::base::Vec3;

use crate::base::{Mole, Stm};
// use crate::orbtprop::obtmat::Mocv;
use numpy::{IntoPyArray, PyArray2, PyArrayMethods};

#[pyclass]
pub struct Calculator {
    pub mole: Mole,
}

impl Calculator {
    // 提取公共的 calculator 创建逻辑
    fn caler(&self) -> rswfn::bondprop::order::Calculator<'_> {
        rswfn::bondprop::order::Calculator::new(&self.mole.core)
    }
}

#[pymethods]
impl Calculator {
    #[new]
    pub fn new(mole: Mole) -> Self {
        Self { mole }
    }

    pub fn mayer(&self, py: Python) -> Py<PyArray2<f64>> {
        let omat = self.caler().mayer();
        let omat = omat.into_pyarray(py).unbind();
        omat
    }

    pub fn wiberg(&self, py: Python) -> Py<PyArray2<f64>> {
        let omat = self.caler().wiberg();
        let omat = omat.into_pyarray(py).unbind();
        omat
    }

    pub fn pocv(
        &self,
        py: Python,
        dirs: HashMap<usize, [f64; 3]>,
        keep_other_atm: bool,
        keep_other_sym: bool,
        otype: &str,
    ) -> Py<PyArray2<f64>> {
        let dirs = dirs
            .into_iter()
            .map(|(key, dir)| (key, Vec3::from_arr(dir)))
            .collect();
        let omat = self
            .caler()
            .pocv(&dirs, keep_other_atm, keep_other_sym, otype);
        omat.into_pyarray(py).unbind()
    }

    pub fn deco(
        &self,
        py: Python,
        stms: HashMap<usize, Stm>, // 每个原子的保留信息
        atos: Vec<usize>,
        otype: &str,
    ) -> Py<PyArray2<f64>> {
        let stms = stms.into_iter().map(|(key, stm)| (key, stm.core)).collect();
        let omat = self.caler().mocv(&stms, &atos, otype);
        omat.into_pyarray(py).unbind()
    }

    pub fn pi_pocv(
        &self,
        py: Python,
        otype: &str,
    ) -> (HashMap<usize, [f64; 3]>, Py<PyArray2<f64>>) {
        let (dirs, omat) = self.caler().pi_pocv(otype);
        let dirs = dirs
            .into_iter()
            .map(|(key, dir)| (key, dir.into_arr()))
            .collect();
        let omat = omat.into_pyarray(py).unbind();
        (dirs, omat)
    }

    pub fn pi_mocv(
        &self,
        py: Python,
        otype: &str,
    ) -> (HashMap<usize, Stm>, Vec<usize>, Py<PyArray2<f64>>) {
        let (stms, atos, omat) = self.caler().pi_mocv(otype);
        let stms = stms
            .into_iter()
            .map(|(key, stm)| (key, Stm { core: stm }))
            .collect();
        let omat = omat.into_pyarray(py).unbind();
        (stms, atos, omat)
    }

    pub fn bound(&self, atm: usize, dir: [f64; 3], nebs: Vec<usize>, otype: &str) -> Vec<f64> {
        let dir = Vec3::from_arr(dir);
        self.caler().bound(atm, &dir, &nebs, otype)
    }
}

#[pyfunction]
pub fn mayer(
    py: Python,
    dmat: Bound<'_, PyArray2<f64>>,
    smat: Bound<'_, PyArray2<f64>>,
) -> Py<PyArray2<f64>> {
    let dmat = dmat.readonly().as_array().to_owned();
    let smat = smat.readonly().as_array().to_owned();
    rswfn::bondprop::order::mayer(&dmat, &smat)
        .into_pyarray(py)
        .unbind()
}

#[pyfunction]
pub fn wiberg(
    py: Python,
    dmat: Bound<'_, PyArray2<f64>>,
    smat: Bound<'_, PyArray2<f64>>,
) -> Py<PyArray2<f64>> {
    let dmat = dmat.readonly().as_array().to_owned();
    let smat = smat.readonly().as_array().to_owned();
    rswfn::bondprop::order::wiberg(&dmat, &smat)
        .into_pyarray(py)
        .unbind()
}

#[pyfunction]
pub fn ato2atm(
    py: Python,
    omat: Bound<'_, PyArray2<f64>>,
    uplo: Vec<[usize; 2]>,
) -> Py<PyArray2<f64>> {
    let omat = omat.readonly().as_array().to_owned();
    rswfn::bondprop::order::ato2atm(&omat, &uplo)
        .into_pyarray(py)
        .unbind()
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "order")?;
    m.add_class::<Calculator>()?;
    m.add_function(wrap_pyfunction!(mayer, &m)?)?;
    m.add_function(wrap_pyfunction!(wiberg, &m)?)?;
    m.add_function(wrap_pyfunction!(ato2atm, &m)?)?;
    parent_module.add_submodule(&m)
}

#[pymodule]
pub mod order {
    #[pymodule_export]
    use super::Calculator;
}
