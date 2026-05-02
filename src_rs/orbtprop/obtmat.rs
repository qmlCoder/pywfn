use pyo3::prelude::*;
use std::collections::HashMap;

use rswfn;
use rswfn::base::Vec3;

use numpy::{IntoPyArray, PyArray2, PyArrayMethods, PyReadonlyArray2};

use crate::base::{Mole, Stm};

#[pyclass]
pub struct Calculator {
    pub mole: Mole,
    pub atms: Vec<usize>,
}

impl Calculator {
    // 提取公共的 calculator 创建逻辑
    fn caler(&self) -> rswfn::orbtprop::obtmat::Calculator<'_> {
        let mut caler = rswfn::orbtprop::obtmat::Calculator::new(&self.mole.core);
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

    pub fn pocv(
        &self,
        py: Python,
        dirs: HashMap<usize, [f64; 3]>, // 每个原子对应的方向向量
        keep_other_atm: bool,           // 其它原子轨道系数是否保留
        keep_other_sym: bool,           //是否保留其它类型的轨道
    ) -> Py<PyArray2<f64>> {
        let dirs = dirs
            .into_iter()
            .map(|(key, dir)| (key, Vec3::from_arr(dir)))
            .collect();
        let cmat_pocv = self.caler().pocv(&dirs, keep_other_atm, keep_other_sym);
        cmat_pocv.into_pyarray(py).unbind()
    }

    pub fn pi_pocv(&self, py: Python) -> (HashMap<usize, [f64; 3]>, Py<PyArray2<f64>>) {
        let (dirs, cmat) = self.caler().pi_pocv();
        let dirs = dirs
            .into_iter()
            .map(|(key, dir)| (key, dir.into_arr()))
            .collect();
        let cmat = cmat.into_pyarray(py).unbind();
        (dirs, cmat)
    }

    pub fn mocv(
        &self,
        py: Python,
        mocvs: HashMap<usize, Mocv>,
        atos: Option<Vec<usize>>,
    ) -> Py<PyArray2<f64>> {
        let mocvs = mocvs
            .into_iter()
            .map(|(key, deco)| (key, deco.core))
            .collect();
        let cmat_deco = self.caler().mocv(&mocvs, atos.as_ref());
        cmat_deco.into_pyarray(py).unbind()
    }

    pub fn pi_mocv(&self, py: Python) -> (HashMap<usize, Mocv>, Py<PyArray2<f64>>) {
        let (mocvs, cmat) = self.caler().pi_mocv();
        let mocvs = mocvs
            .into_iter()
            .map(|(key, mocv)| (key, Mocv { core: mocv }))
            .collect();
        let cmat = cmat.into_pyarray(py).unbind();
        (mocvs, cmat)
    }

    /// 局部坐标系的系数矩阵
    pub fn ldao(&self, py: Python, stms: HashMap<usize, Stm>) -> Py<PyArray2<f64>> {
        let stms = stms.into_iter().map(|(idx, stm)| (idx, stm.core)).collect();
        let cmat = self.caler().ldao(&stms);
        let cmat = cmat.into_pyarray(py).unbind();
        cmat
    }
}

#[pyclass]
#[derive(Clone)]
pub struct Mocv {
    pub core: rswfn::orbtprop::obtmat::Mocv,
}

#[pymethods]
impl Mocv {
    #[new]
    pub fn new(
        _py: Python,
        stm: Stm,
        skeep: [usize; 1],
        pkeep: [usize; 3],
        dkeep: [usize; 5],
    ) -> Self {
        let stm = stm.core;
        Self {
            core: rswfn::orbtprop::obtmat::Mocv {
                stm,
                skeep,
                pkeep,
                dkeep,
            },
        }
    }

    pub fn __repr__(&self) -> String {
        format!("{}", self.core)
    }
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "obtmat")?;
    m.add_class::<Calculator>()?;
    m.add_class::<Mocv>()?;
    parent_module.add_submodule(&m)
}

#[pymodule]
pub mod obtmat {
    #[pymodule_export]
    use super::Calculator;
    #[pymodule_export]
    use super::Mocv;
}
