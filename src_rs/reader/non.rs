#![allow(non_snake_case)]

use pyo3::prelude::*;
use rswfn;

use crate::base::Atoms;
use rswfn::reader::Reader;

#[pyclass]
pub struct NonReader {
    core: rswfn::reader::NonReader,
}

#[pymethods]
impl NonReader {
    #[new]
    pub fn new(path: String) -> PyResult<Self> {
        let result = Self {
            core: rswfn::reader::NonReader::new(&path),
        };
        Ok(result)
    }

    pub fn set_atoms(&mut self, atms: Vec<usize>, xyzs: Vec<[f64; 3]>) {
        self.core.set_atoms(&atms, &xyzs);
    }
}
