#![allow(non_snake_case)]

use numpy::{IntoPyArray, PyArray2};
use pyo3::prelude::*;
use rswfn;

use crate::base::{Atoms, Basis, Coefs};
use crate::reader::Nics;
use rswfn::reader::Reader;

#[pyclass]
pub struct LogReader {
    core: rswfn::reader::LogReader,
}

#[pymethods]
impl LogReader {
    #[new]
    pub fn new(path: String) -> PyResult<Self> {
        let reader = Self {
            core: rswfn::reader::LogReader::new(&path),
        };
        Ok(reader)
    }

    pub fn read_job(&self) -> String {
        self.core.read_job().to_owned()
    }

    pub fn read_smat<'py>(&self, py: Python<'py>) -> Option<Bound<'py, PyArray2<f64>>> {
        self.core.read_smat().map(|smat| smat.into_pyarray(py))
    }

    pub fn read_tmat<'py>(&self, py: Python<'py>) -> Option<Bound<'py, PyArray2<f64>>> {
        self.core.read_tmat().map(|tmat| tmat.into_pyarray(py))
    }

    pub fn read_vmat<'py>(&self, py: Python<'py>) -> Option<Bound<'py, PyArray2<f64>>> {
        self.core.read_vmat().map(|vmat| vmat.into_pyarray(py))
    }

    pub fn read_nics<'py>(&self, py: Python<'py>) -> Option<Vec<Nics>> {
        if let Some(nics) = self.core.read_nics() {
            Some(nics.into_iter().map(|core| Nics { core }).collect())
        } else {
            None
        }
    }

    pub fn get_atoms(&self) -> Atoms {
        let atoms = self.core.get_atoms();
        Atoms { core: atoms }
    }

    pub fn get_basis(&self) -> Option<Basis> {
        let basis = self.core.get_basis();
        if let Some(basis) = basis {
            Some(Basis { core: basis })
        } else {
            None
        }
    }

    pub fn get_coefs(&self) -> Option<Coefs> {
        let coefs = self.core.get_coefs();
        if let Some(coefs) = coefs {
            Some(Coefs { core: coefs })
        } else {
            None
        }
    }
}
