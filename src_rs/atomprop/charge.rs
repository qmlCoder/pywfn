use ndarray::{Array, Array1, Array2};
use pyo3::prelude::*;
use std::collections::HashMap;

use rswfn;
use rswfn::base::Vec3;

use crate::base::{Mole, Stm};
// use crate::orbtprop::obtmat::Mocv;
use numpy::{IntoPyArray, PyArray1, PyArray2, PyArrayMethods};

#[pyclass]
pub struct Calculator {
    pub mole: Mole,
    pub atms: Vec<usize>,
}

impl Calculator {
    // 提取公共的 calculator 创建逻辑
    fn caler(&self) -> rswfn::atomprop::charge::Calculator<'_> {
        let mut caler = rswfn::atomprop::charge::Calculator::new(&self.mole.core);
        caler.atms = &self.atms;
        caler
    }
}

#[pymethods]
impl Calculator {
    #[new]
    pub fn new(mole: Mole, atms: Option<Vec<usize>>) -> Self {
        let atms = match atms {
            Some(atms) => atms,
            None => mole.core.atoms().idxs().clone(),
        };
        Self { mole, atms }
    }

    pub fn spin(&self, py: Python, ctype: &str) -> Py<PyArray2<f64>> {
        let vals = self.caler().spin(ctype);
        let arr = Array2::from_shape_vec((vals.len(), 3), vals.concat()).unwrap();
        arr.into_pyarray(py).unbind()
    }

    pub fn mulliken(&self, py: Python) -> Py<PyArray1<f64>> {
        let vals = self.caler().mulliken();
        let arr = Array1::from_vec(vals);
        arr.into_pyarray(py).unbind()
    }

    pub fn lowdin(&self, py: Python) -> Py<PyArray1<f64>> {
        let vals = self.caler().lowdin();
        let arr = Array1::from_vec(vals);
        arr.into_pyarray(py).unbind()
    }

    pub fn hirshfeld(&self, py: Python) -> Py<PyArray1<f64>> {
        let vals = self.caler().hirshfeld();
        let arr = Array1::from_vec(vals);
        arr.into_pyarray(py).unbind()
    }

    pub fn pocv(
        &self,
        dirs: HashMap<usize, [f64; 3]>,
        keep_other_atm: bool,
        keep_other_sym: bool,
        ctype: &str,
    ) -> PyResult<Vec<f64>> {
        let dirs: HashMap<usize, Vec3> = dirs
            .into_iter()
            .map(|(atm, dir)| (atm, Vec3::from_arr(dir)))
            .collect();
        let vals = self
            .caler()
            .pocv(&dirs, keep_other_atm, keep_other_sym, ctype);
        Ok(vals)
    }

    pub fn mocv(
        &self,
        stms: HashMap<usize, Stm>,
        atos: Vec<usize>,
        ctype: &str,
    ) -> PyResult<Vec<f64>> {
        let stms = stms
            .into_iter()
            .map(|(atm, stm)| (atm, stm.core.clone()))
            .collect();
        let vals = self.caler().mocv(&stms, &atos, ctype);
        Ok(vals)
    }

    pub fn pi_pocv(&self, ctype: &str) -> PyResult<(HashMap<usize, [f64; 3]>, Vec<f64>)> {
        let (dirs, eles) = self.caler().pi_pocv(ctype);
        let dirs = dirs
            .iter()
            .map(|(key, dir)| (*key, dir.into_arr()))
            .collect();
        Ok((dirs, eles))
    }

    pub fn pi_mocv(&self, ctype: &str) -> PyResult<(HashMap<usize, Stm>, Vec<usize>, Vec<f64>)> {
        let (stms, atos, eles) = self.caler().pi_mocv(ctype);
        let stms = stms
            .iter()
            .map(|(atm, stm)| (*atm, Stm { core: stm.clone() }))
            .collect();
        Ok((stms, atos, eles))
    }
}

#[pyfunction]
pub fn mulliken(
    dmat: Bound<'_, PyArray2<f64>>,
    smat: Bound<'_, PyArray2<f64>>,
    uplo: Vec<[usize; 2]>,
) -> Vec<f64> {
    let dmat = dmat.readonly().as_array().to_owned();
    let smat = smat.readonly().as_array().to_owned();
    rswfn::atomprop::charge::mulliken(&dmat, &smat, &uplo)
}

#[pyfunction]
pub fn lowdin(
    dmat: Bound<'_, PyArray2<f64>>,
    smat: Bound<'_, PyArray2<f64>>,
    uplo: Vec<[usize; 2]>,
) -> Vec<f64> {
    let dmat = dmat.readonly().as_array().to_owned();
    let smat = smat.readonly().as_array().to_owned();
    rswfn::atomprop::charge::lowdin(&dmat, &smat, &uplo)
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "charge")?;
    m.add_class::<Calculator>()?;
    m.add_function(wrap_pyfunction!(mulliken, &m)?)?;
    m.add_function(wrap_pyfunction!(lowdin, &m)?)?;
    parent_module.add_submodule(&m)
}

#[pymodule]
pub mod charge {
    #[pymodule_export]
    use super::Calculator;
}
