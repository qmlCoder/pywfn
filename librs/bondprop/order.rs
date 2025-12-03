use nalgebra::Vector3;
use pyo3::prelude::*;
use std::collections::HashMap;

use rswfn;

use crate::base::Mole;
use crate::orbtprop::obtmat::Deco;
use numpy::{IntoPyArray, PyArray2};

#[pyclass]
pub struct Calculator {
    pub mole: Mole,
}

impl Calculator {
    // 提取公共的 calculator 创建逻辑
    fn caler(&self) -> rswfn::bondprop::order::Calculator<'_> {
        rswfn::bondprop::order::Calculator::new(&self.mole.inner)
    }
}

#[pymethods]
impl Calculator {
    #[new]
    pub fn new(mole: Mole) -> Self {
        Self { mole }
    }

    pub fn mayer(&self, py: Python) -> Py<PyArray2<f64>> {
        let orders = self.caler().mayer();
        let orders = orders.into_pyarray(py).unbind();
        orders
    }

    pub fn pocv(
        &self,
        py: Python,
        dirs: HashMap<u32, [f64; 3]>,
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
        decos: HashMap<u32, Deco>, // 每个原子的保留信息
        otype: &str,
    ) -> Py<PyArray2<f64>> {
        let decos = decos
            .into_iter()
            .map(|(key, deco)| (key, deco.inner))
            .collect();
        let orders = self.caler().deco(&decos, otype);
        orders.into_pyarray(py).unbind()
    }

    pub fn pi_pocv(&self, py: Python) -> (HashMap<u32, [f64; 3]>, Py<PyArray2<f64>>) {
        let (dirs, omat) = self.caler().pi_pocv();
        // let nrow = omat.nrows();
        // let ncol = omat.ncols();
        // for i in 0..nrow {
        //     for j in 0..ncol {
        //         if i == j {
        //             continue;
        //         }
        //         if omat[(i, j)] < 1e-3 {
        //             continue;
        //         }
        //         println!("{:>2}-{:>2}:{:>10.3}", i + 1, j + 1, omat[(i, j)]);
        //     }
        // }
        let dirs = dirs
            .into_iter()
            .map(|(key, dir)| (key, [dir.x, dir.y, dir.z]))
            .collect();
        let omat = omat.into_pyarray(py).unbind();
        (dirs, omat)
    }

    pub fn pi_deco(&self, py: Python) -> (HashMap<u32, Deco>, Py<PyArray2<f64>>) {
        let (decos, omat) = self.caler().pi_deco();
        let decos = decos
            .into_iter()
            .map(|(key, deco)| (key, Deco { inner: deco }))
            .collect();
        let omat = omat.into_pyarray(py).unbind();
        (decos, omat)
    }

    pub fn bound(&self, atm: u32, dir: [f64; 3], otype: &str) -> (Vec<u32>, Vec<f64>) {
        let dir = Vector3::new(dir[0], dir[1], dir[2]);
        self.caler().bound(atm, &dir, otype)
    }
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "order")?;
    m.add_class::<Calculator>()?;
    parent_module.add_submodule(&m)
}
