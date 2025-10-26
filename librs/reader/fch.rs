#![allow(non_snake_case)]

use pyo3::prelude::*;
use rswfn;

use crate::base::{Basis, Coefs, Geome};
use rswfn::reader::Reader;

#[pyclass]
pub struct FchReader {
    inner: rswfn::reader::FchReader,
}

#[pymethods]
impl FchReader {
    #[new]
    pub fn new(path: String) -> PyResult<Self> {
        let reader = Self {
            inner: rswfn::reader::FchReader::new(&path),
        };
        Ok(reader)
    }

    pub fn get_geome(&self) -> Geome {
        let geome = self.inner.get_geome();
        Geome { inner: geome }
    }

    pub fn get_basis(&self) -> Basis {
        let basis = self.inner.get_basis();
        Basis { inner: basis }
    }

    pub fn get_coefs(&self) -> Coefs {
        let coefs = self.inner.get_coefs();
        Coefs { inner: coefs }
    }
}
