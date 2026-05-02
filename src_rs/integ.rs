#![allow(non_snake_case)]

use std::collections::HashMap;

use crate::base::{Mole, Stm};
use numpy::{IntoPyArray, PyArray2, PyArray4, PyReadonlyArray2};
use pyo3::prelude::*;
use rswfn;

#[pyclass]
pub struct Int_mat {
    pub mole: Mole,
}

impl Int_mat {
    // 提取公共的 calculator 创建逻辑
    fn caler(&self) -> rswfn::integ::Int_mat<'_> {
        rswfn::integ::Int_mat::new(&self.mole.core)
    }
}

#[pymethods]
impl Int_mat {
    #[new]
    pub fn new(mole: Mole) -> Self {
        Self { mole }
    }

    pub fn smat(&self, py: Python, form: String) -> Py<PyArray2<f64>> {
        let smat = self.caler().smat(&form).clone();
        smat.into_pyarray(py).unbind()
    }

    pub fn tmat(&self, py: Python, form: String) -> Py<PyArray2<f64>> {
        let tmat = self.caler().tmat(&form).clone();
        tmat.into_pyarray(py).unbind()
    }

    pub fn vmat(&self, py: Python, form: String) -> Py<PyArray2<f64>> {
        let vmat = self.caler().vmat(&form).clone();
        vmat.into_pyarray(py).unbind()
    }

    pub fn mat1e(&self, py: Python, cals: [usize; 3], form: String) -> [Py<PyArray2<f64>>; 3] {
        let [smat, tmat, vmat] = self.caler().mat1e(cals, &form);
        let smat = smat.into_pyarray(py).unbind();
        let tmat = tmat.into_pyarray(py).unbind();
        let vmat = vmat.into_pyarray(py).unbind();
        [smat, tmat, vmat]
    }

    pub fn mat2e(&self, py: Python, form: String) -> Py<PyArray4<f64>> {
        let mat2e = self.caler().mat2e(&form);
        let mat2e = mat2e.into_mat4();
        mat2e.into_pyarray(py).unbind()
    }

    pub fn mat1e_lcs(
        &self,
        py: Python,
        stms: HashMap<usize, Stm>,
        cals: [usize; 3],
    ) -> [Py<PyArray2<f64>>; 3] {
        let stms = stms.into_iter().map(|(idx, stm)| (idx, stm.core)).collect();
        let [smat, tmat, vmat] = self.caler().mat1e_lcs(&stms, cals);
        let smat = smat.into_pyarray(py).unbind();
        let tmat = tmat.into_pyarray(py).unbind();
        let vmat = vmat.into_pyarray(py).unbind();
        [smat, tmat, vmat]
    }
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "integ")?; // 为当前的rs文件创建一个子模块
    m.add_class::<Int_mat>()?;
    // m.add_function(wrap_pyfunction!(rys_roots, &m)?)?;
    parent_module.add_submodule(&m) // 将子模块添加到父模块中
}
