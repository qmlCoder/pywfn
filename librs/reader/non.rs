#![allow(non_snake_case)]

use pyo3::prelude::*;
use rswfn;

use crate::base::Geome;
use rswfn::reader::Reader;

#[pyclass]
pub struct NonReader {
    inner: rswfn::reader::NonReader,
}

#[pymethods]
impl NonReader {
    #[new]
    pub fn new(path: String) -> PyResult<Self> {
        let result = Self {
            inner: rswfn::reader::NonReader::new(&path),
        };
        Ok(result)
    }

    pub fn set_geome(&mut self, atms: Vec<u32>, xyzs: Vec<[f64; 3]>) {
        self.inner.set_geome(&atms, &xyzs);
    }
}
