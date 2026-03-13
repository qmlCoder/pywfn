#![allow(non_snake_case)]

use pyo3::prelude::*;
use rswfn;

use crate::base::{Basis, Coefs, Geome};
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

    pub fn get_geome(&self) -> Geome {
        let geome = self.core.get_geome();
        Geome { core: geome }
    }

    pub fn get_basis(&self) -> Basis {
        let basis = self.core.get_basis();
        Basis { core: basis }
    }

    pub fn get_coefs(&self) -> Coefs {
        let coefs = self.core.get_coefs();
        Coefs { core: coefs }
    }

    pub fn get_neles(&self) -> [usize; 2] {
        self.core.get_neles()
    }
}
