#![allow(non_snake_case)]

use pyo3::prelude::*;
use rswfn;

use crate::base::{Atoms, Basis, Coefs};
use rswfn::reader::Reader;

#[pyclass]
pub struct CubReader {
    core: rswfn::reader::CubReader,
}

#[pymethods]
impl CubReader {
    #[new]
    pub fn new(path: String) -> PyResult<Self> {
        let reader = Self {
            core: rswfn::reader::CubReader::new(&path),
        };
        Ok(reader)
    }

    pub fn get_atoms(&self) -> Atoms {
        let atoms = self.core.get_atoms();
        Atoms { core: atoms }
    }

    pub fn get_basis(&self) -> Basis {
        let basis = self.core.get_basis();
        Basis { core: basis }
    }

    pub fn get_coefs(&self) -> Coefs {
        let coefs = self.core.get_coefs();
        Coefs { core: coefs }
    }
}
