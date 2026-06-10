use numpy::{IntoPyArray, PyArray2};
use pyo3::prelude::*;

use rswfn;

use crate::base::Atoms;

#[pyclass]
pub struct Atos {
    pub core: rswfn::base::Atos,
}

#[pymethods]
impl Atos {
    pub fn len(&self, form: &str) -> usize {
        self.core.len(form)
    }

    pub fn atms(&self, form: &str) -> Vec<usize> {
        self.core.atms(form)
    }

    pub fn shls(&self, form: &str) -> Vec<usize> {
        self.core.shls(form)
    }

    pub fn syms(&self, form: &str) -> Vec<String> {
        self.core.syms(form)
    }
}

#[pyclass]
pub struct Shl {
    pub core: rswfn::base::Shl,
}

#[pymethods]
impl Shl {
    #[getter]
    pub fn atm(&self) -> usize {
        self.core.atm
    }
    #[getter]
    pub fn shl(&self) -> usize {
        self.core.shl
    }
    #[getter]
    pub fn ang(&self) -> usize {
        self.core.ang
    }
    #[getter]
    pub fn loc(&self) -> usize {
        self.core.loc()
    }
    #[getter]
    pub fn form(&self) -> String {
        self.core.form.clone()
    }
    #[getter]
    pub fn alps(&self) -> Vec<f64> {
        self.core.alps.clone()
    }
    #[getter]
    pub fn coes(&self) -> Vec<f64> {
        self.core.coes.clone()
    }
}

#[pyclass]
pub struct Shls {
    pub core: rswfn::base::Shls,
}

#[pymethods]
impl Shls {
    pub fn get(&self, i: usize) -> Shl {
        Shl {
            core: self.core.get(i).clone(),
        }
    }

    pub fn len(&self) -> usize {
        self.core.len()
    }

    pub fn uplo(&self, form: &str) -> Vec<[usize; 2]> {
        self.core.uplo(form).clone()
    }
}

#[derive(Clone)]
#[pyclass]
pub struct Basis {
    pub core: rswfn::base::Basis,
}

#[pymethods]
impl Basis {
    pub fn __repr__(&self) -> String {
        format!("{}", self.core)
    }

    pub fn atos(&self) -> Atos {
        Atos {
            core: self.core.atos().clone(),
        }
    }

    pub fn shls(&self) -> Shls {
        Shls {
            core: self.core.shls.clone(),
        }
    }

    pub fn get_smat(&self, py: Python, form: &str, atoms: Atoms) -> Py<PyArray2<f64>> {
        let smat = self.core.get_smat(form, &atoms.core).to_owned();
        smat.into_pyarray(py).unbind()
    }

    pub fn atm_uplo(&self, form: &str) -> Vec<[usize; 2]> {
        let uplo = self.core.atm_uplo(form).to_owned();
        uplo
    }
}
