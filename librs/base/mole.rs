use numpy::{IntoPyArray, PyArray2};
use pyo3::prelude::*;

use crate::base::{Basis, Coefs, Geome};
use rswfn;

#[derive(Clone)]
#[pyclass]
pub struct Mole {
    pub inner: rswfn::base::Mole,
}

#[pymethods]
impl Mole {
    #[new]
    pub fn new(geome: Geome, basis: Basis, coefs: Coefs) -> PyResult<Mole> {
        let geome = Some(geome.inner);
        let basis = Some(basis.inner);
        let coefs = Some(coefs.inner);
        let mole = Mole {
            inner: rswfn::base::Mole::new(geome, basis, coefs),
        };
        Ok(mole)
    }

    pub fn geome(&self) -> Geome {
        let geome = self.inner.geome().clone();
        Geome { inner: geome }
    }

    pub fn syms(&self) -> Vec<String> {
        self.inner.atoms().syms()
    }

    pub fn xyzs(&self, angstrom: bool) -> Vec<[f64; 3]> {
        self.inner.atoms().xyzs(angstrom)
    }

    pub fn get_cmat(&self, py: Python, form: &str) -> Py<PyArray2<f64>> {
        let cmat = self.inner.coefs().get_cmat(form).clone();
        cmat.into_pyarray(py).unbind()
    }

    pub fn get_smat(&self, py: Python, form: &str) -> Py<PyArray2<f64>> {
        let smat = self.inner.coefs().get_smat(form, &self.inner).clone();
        smat.into_pyarray(py).unbind()
    }

    pub fn get_dmat(&self, py: Python, form: &str) -> Py<PyArray2<f64>> {
        let dmat = self.inner.coefs().get_dmat(form).clone();
        dmat.into_pyarray(py).unbind()
    }

    pub fn is_open(&self) -> bool {
        self.inner.is_open()
    }
}
