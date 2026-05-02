use pyo3::prelude::*;
use std::collections::HashMap;

use rswfn;
use rswfn::base::Vec3;

use crate::base::Mole;
use crate::orbtprop::obtmat::Mocv;
use numpy::{IntoPyArray, PyArray2};

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

    pub fn lowdin(&self, py: Python) -> Py<PyArray2<f64>> {
        let omat = self.caler().lowdin();
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
        mocvs: HashMap<usize, Mocv>, // 每个原子的保留信息
        otype: &str,
    ) -> Py<PyArray2<f64>> {
        let mocvs = mocvs
            .into_iter()
            .map(|(key, mocv)| (key, mocv.core))
            .collect();
        let omat = self.caler().mocv(&mocvs, otype);
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

    pub fn pi_mocv(&self, py: Python, otype: &str) -> (HashMap<usize, Mocv>, Py<PyArray2<f64>>) {
        let (mocvs, omat) = self.caler().pi_mocv(otype);
        let mocvs = mocvs
            .into_iter()
            .map(|(key, mocv)| (key, Mocv { core: mocv }))
            .collect();
        let omat = omat.into_pyarray(py).unbind();
        (mocvs, omat)
    }

    pub fn bound(&self, atm: usize, dir: [f64; 3], nebs: Vec<usize>, otype: &str) -> Vec<f64> {
        let dir = Vec3::from_arr(dir);
        self.caler().bound(atm, &dir, &nebs, otype)
    }
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "order")?;
    m.add_class::<Calculator>()?;
    parent_module.add_submodule(&m)
}

#[pymodule]
pub mod order {
    #[pymodule_export]
    use super::Calculator;
}
