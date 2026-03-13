use numpy::{IntoPyArray, PyArray2};
use pyo3::prelude::*;

use crate::base::{Basis, Coefs, Geome};
use rswfn;

#[derive(Clone)]
#[pyclass]
pub struct Mole {
    pub core: rswfn::base::Mole,
}

#[pymethods]
impl Mole {
    #[new]
    pub fn new(geome: Geome, basis: Basis, coefs: Coefs, neles: [usize; 2]) -> PyResult<Mole> {
        let geome = geome.core;
        let basis = basis.core;
        let coefs = coefs.core;
        let mole = Mole {
            core: rswfn::base::Mole::new(geome, basis, coefs, neles),
        };
        Ok(mole)
    }

    pub fn geome(&self) -> Geome {
        let geome = self.core.geome().clone();
        Geome { core: geome }
    }

    pub fn syms(&self) -> Vec<String> {
        self.core.atoms().syms()
    }

    pub fn xyzs(&self) -> Vec<[f64; 3]> {
        self.core.atoms().xyzs().clone()
    }

    pub fn get_cmat(&self, py: Python, form: &str) -> Py<PyArray2<f64>> {
        let cmat = self.core.coefs().get_cmat(form).clone();
        cmat.into_pyarray(py).unbind()
    }

    pub fn get_smat(&self, py: Python, form: &str) -> Py<PyArray2<f64>> {
        let smat = self.core.coefs().get_smat(form, &self.core).clone();
        smat.into_pyarray(py).unbind()
    }

    pub fn get_dmat(&self, py: Python, form: &str) -> Py<PyArray2<f64>> {
        let dmat = self.core.coefs().get_dmat(form).clone();
        dmat.into_pyarray(py).unbind()
    }

    pub fn is_open(&self) -> bool {
        self.core.is_open()
    }
}
