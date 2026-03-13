use nalgebra::Vector3;
use pyo3::prelude::*;
use std::collections::HashMap;

use rswfn;

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
        let orders = self.caler().order("mayer");
        let orders = orders.into_pyarray(py).unbind();
        orders
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
            .map(|(key, dir)| (key, Vector3::new(dir[0], dir[1], dir[2])))
            .collect();
        let orders = self
            .caler()
            .pocv(&dirs, keep_other_atm, keep_other_sym, otype);
        orders.into_pyarray(py).unbind()
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
        let orders = self.caler().mocv(&mocvs, otype);
        orders.into_pyarray(py).unbind()
    }

    pub fn pi_pocv(&self, py: Python) -> (HashMap<usize, [f64; 3]>, Py<PyArray2<f64>>) {
        let (dirs, omat) = self.caler().pi_pocv();
        let dirs = dirs
            .into_iter()
            .map(|(key, dir)| (key, [dir.x, dir.y, dir.z]))
            .collect();
        let omat = omat.into_pyarray(py).unbind();
        (dirs, omat)
    }

    pub fn pi_mocv(&self, py: Python) -> (HashMap<usize, Mocv>, Py<PyArray2<f64>>) {
        let (mocvs, omat) = self.caler().pi_mocv();
        let mocvs = mocvs
            .into_iter()
            .map(|(key, mocv)| (key, Mocv { core: mocv }))
            .collect();
        let omat = omat.into_pyarray(py).unbind();
        (mocvs, omat)
    }

    pub fn bound(&self, atm: usize, dir: [f64; 3], otype: &str) -> (Vec<usize>, Vec<f64>) {
        let dir = Vector3::new(dir[0], dir[1], dir[2]);
        self.caler().bound(atm, &dir, otype)
    }
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "order")?;
    m.add_class::<Calculator>()?;
    parent_module.add_submodule(&m)
}
