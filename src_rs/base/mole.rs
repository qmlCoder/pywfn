use crate::base::{Atoms, Basis, Coefs};
use numpy::PyArrayMethods;
use numpy::{IntoPyArray, PyArray2};
use pyo3::prelude::*;
use rswfn;

#[derive(Clone)]
#[pyclass]
pub struct Mole {
    pub core: rswfn::base::Mole,
}

#[pymethods]
impl Mole {
    #[new]
    pub fn new(path: String) -> PyResult<Mole> {
        let mole = Mole {
            core: rswfn::base::Mole::from_file(&path),
        };
        Ok(mole)
    }

    pub fn atoms(&self) -> Atoms {
        let atoms = self.core.atoms().clone();
        Atoms { core: atoms }
    }

    pub fn basis(&self) -> Basis {
        let basis = self.core.basis().clone();
        Basis { core: basis }
    }

    pub fn coefs(&self) -> Coefs {
        let coefs = self.core.coefs().clone();
        Coefs { core: coefs }
    }

    pub fn syms(&self) -> Vec<String> {
        self.core.atoms().syms()
    }

    pub fn xyzs(&self) -> Vec<[f64; 3]> {
        self.core.atoms().xyzs().clone()
    }

    pub fn set_cmat(&mut self, form: &str, cmat: Bound<'_, PyArray2<f64>>) {
        let cmat = cmat.readonly().as_array().to_owned();
        self.core.coefs.set_cmat(form, &cmat);
    }

    pub fn get_cmat(&self, py: Python, form: &str) -> Py<PyArray2<f64>> {
        let cmat = self.core.coefs().get_cmat(form).clone();
        cmat.into_pyarray(py).unbind()
    }

    pub fn get_smat(&self, py: Python, form: &str) -> Py<PyArray2<f64>> {
        let smat = self.core.basis().get_smat(form, &self.core).clone();
        smat.into_pyarray(py).unbind()
    }

    pub fn get_dmat(&self, py: Python, form: &str) -> Py<PyArray2<f64>> {
        let dmat = self.core.coefs().get_dmat(form).clone();
        dmat.into_pyarray(py).unbind()
    }

    pub fn is_open(&self) -> bool {
        self.core.is_open()
    }

    pub fn border(&self) -> [[f64; 3]; 2] {
        self.core.border()
    }

    pub fn __repr__(&self) -> String {
        format!("{}", self.core)
    }
}
