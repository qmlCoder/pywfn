use nalgebra::Vector3;
use pyo3::prelude::*;
use std::collections::HashMap;

use rswfn;

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
        let mut caler = rswfn::orbtprop::obtmat::Calculator::new(&self.mole.inner);
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
            None => mole.inner.atoms().idxs().clone(),
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
            .map(|(key, dir)| (key, Vector3::new(dir[0], dir[1], dir[2])))
            .collect();
        let cmat_pocv = self.caler().pocv(&dirs, keep_other_atm, keep_other_sym);
        cmat_pocv.into_pyarray(py).unbind()
    }

    pub fn pi_pocv(&self, py: Python) -> (HashMap<usize, [f64; 3]>, Py<PyArray2<f64>>) {
        let (dirs, cmat) = self.caler().pi_pocv();
        let dirs = dirs
            .into_iter()
            .map(|(key, dir)| (key, [dir.x, dir.y, dir.z]))
            .collect();
        let cmat = cmat.into_pyarray(py).unbind();
        (dirs, cmat)
    }

    pub fn deco(&self, py: Python, decos: HashMap<usize, Deco>) -> Py<PyArray2<f64>> {
        let decos = decos
            .into_iter()
            .map(|(key, deco)| (key, deco.inner))
            .collect();
        let cmat_deco = self.caler().deco(&decos);
        cmat_deco.into_pyarray(py).unbind()
    }

    pub fn pi_deco(&self, py: Python) -> (HashMap<usize, Deco>, Py<PyArray2<f64>>) {
        let (decos, cmat) = self.caler().pi_deco();
        let decos = decos
            .into_iter()
            .map(|(key, deco)| (key, Deco { inner: deco }))
            .collect();
        let cmat = cmat.into_pyarray(py).unbind();
        (decos, cmat)
    }

    /// 局部坐标系的系数矩阵
    pub fn ldao(&self, py: Python, stms: Vec<Stm>) -> Py<PyArray2<f64>> {
        let stms = stms.into_iter().map(|stm| stm.inner).collect();
        let cmat = self.caler().ldao(&stms);
        let cmat = cmat.into_pyarray(py).unbind();
        cmat
    }
}

#[pyclass]
#[derive(Clone)]
pub struct Deco {
    pub inner: rswfn::orbtprop::obtmat::Deco,
}

#[pymethods]
impl Deco {
    #[new]
    pub fn new(
        _py: Python,
        stm: PyReadonlyArray2<f64>,
        skeep: [usize; 1],
        pkeep: [usize; 3],
        dkeep: [usize; 6],
    ) -> Self {
        let stm = stm.to_owned_array();
        let stm = rswfn::base::Stm::from_mat(stm);
        Self {
            inner: rswfn::orbtprop::obtmat::Deco {
                stm,
                skeep,
                pkeep,
                dkeep,
            },
        }
    }

    pub fn __repr__(&self) -> String {
        format!("{}", self.inner)
    }
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "obtmat")?;
    m.add_class::<Calculator>()?;
    m.add_class::<Deco>()?;
    parent_module.add_submodule(&m)
}
