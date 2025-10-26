use nalgebra::Vector3;
use pyo3::prelude::*;
use std::collections::HashMap;

use rswfn;

use numpy::{IntoPyArray, PyArray2, PyArrayMethods, PyReadonlyArray2};

use crate::base::Mole;

#[pyclass]
pub struct Calculator {
    pub mole: Mole,
}

impl Calculator {
    // 提取公共的 calculator 创建逻辑
    fn caler(&self) -> rswfn::moleprop::orbital::Calculator<'_> {
        rswfn::moleprop::orbital::Calculator::new(&self.mole.inner)
    }
}

#[pymethods]
impl Calculator {
    #[new]
    pub fn new(mole: Mole) -> Self {
        Self { mole }
    }

    pub fn pocv(
        &self,
        py: Python,
        dirs: HashMap<u32, [f64; 3]>, // 每个原子对应的方向向量
        keep_other_atm: bool,         // 其它原子轨道系数是否保留
        keep_other_sym: bool,         //是否保留其它类型的轨道
    ) -> Py<PyArray2<f64>> {
        let dirs = dirs
            .into_iter()
            .map(|(key, dir)| (key, Vector3::new(dir[0], dir[1], dir[2])))
            .collect();
        let cmat_pocv = self.caler().pocv(&dirs, keep_other_atm, keep_other_sym);
        cmat_pocv.into_pyarray(py).unbind()
    }

    pub fn pi_pocv(
        &self,
        py: Python,
        atms: Vec<u32>,
    ) -> (HashMap<u32, [f64; 3]>, Py<PyArray2<f64>>) {
        let atms = if atms.len() == 0 { None } else { Some(&atms) };
        let (dirs, cmat) = self.caler().pi_pocv(atms);
        let dirs = dirs
            .into_iter()
            .map(|(key, dir)| (key, [dir.x, dir.y, dir.z]))
            .collect();
        let cmat = cmat.into_pyarray(py).unbind();
        (dirs, cmat)
    }

    pub fn deco(&self, py: Python, decos: HashMap<u32, Deco>) -> Py<PyArray2<f64>> {
        let decos = decos
            .into_iter()
            .map(|(key, deco)| (key, deco.inner))
            .collect();
        let cmat_deco = self.caler().deco(&decos);
        cmat_deco.into_pyarray(py).unbind()
    }

    pub fn pi_deco(&self, py: Python) -> (HashMap<u32, Deco>, Py<PyArray2<f64>>) {
        let (decos, cmat) = self.caler().pi_deco();
        let decos = decos
            .into_iter()
            .map(|(key, deco)| (key, Deco { inner: deco }))
            .collect();
        let cmat = cmat.into_pyarray(py).unbind();
        (decos, cmat)
    }
}

#[pyclass]
#[derive(Clone)]
pub struct Deco {
    pub inner: rswfn::moleprop::orbital::Deco,
}

#[pymethods]
impl Deco {
    #[new]
    pub fn new(
        _py: Python,
        base: PyReadonlyArray2<f64>,
        skeep: [u32; 1],
        pkeep: [u32; 3],
        dkeep: [u32; 6],
    ) -> Self {
        let base = base.to_owned_array();
        Self {
            inner: rswfn::moleprop::orbital::Deco {
                base,
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
    let m = PyModule::new(parent_module.py(), "orbital")?;
    m.add_class::<Calculator>()?;
    m.add_class::<Deco>()?;
    parent_module.add_submodule(&m)
}
