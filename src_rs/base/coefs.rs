use pyo3::prelude::*;

use numpy::{IntoPyArray, PyArray2, PyReadonlyArray2};
use rswfn;

use crate::{
    base::{Basis, Mole},
    utils,
};

#[derive(Clone)]
#[pyclass]
pub struct Coefs {
    pub core: rswfn::base::Coefs,
}

#[pymethods]
impl Coefs {
    pub fn get_cmat(&self, py: Python, form: &str, basis: Basis) -> Py<PyArray2<f64>> {
        let cmat = self.core.get_cmat(form, &basis.core).to_owned();
        cmat.into_pyarray(py).unbind()
    }

    pub fn get_dmat(&self, py: Python, form: &str, basis: Basis) -> Py<PyArray2<f64>> {
        let dmat = self.core.get_dmat(form, &basis.core).to_owned();
        dmat.into_pyarray(py).unbind()
    }

    #[getter]
    pub fn obt_occs(&self) -> Vec<f64> {
        self.core.obt_occs.clone()
    }

    #[getter]
    pub fn obt_engs(&self) -> Vec<f64> {
        self.core.obt_engs.clone()
    }

    pub fn obt_idxs(&self, if_occ: bool) -> Vec<usize> {
        self.core.obt_idxs(if_occ).clone()
    }

    pub fn __repr__(&self) -> String {
        format!("{}", self.core)
    }
}
