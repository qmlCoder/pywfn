use pyo3::prelude::*;

use rswfn;

#[derive(Clone)]
#[pyclass]
pub struct Geome {
    pub inner: rswfn::base::Geome,
}

#[pymethods]
impl Geome {
    #[new]
    pub fn new() -> PyResult<Self> {
        Ok(Self {
            inner: rswfn::base::Geome::blank(),
        })
    }

    pub fn build(&mut self, atms: Vec<u32>, xyzs: Vec<[f64; 3]>, crgs: Vec<f64>) {
        self.inner.build(&atms, &xyzs, &crgs);
    }

    pub fn __repr__(&self) -> String {
        format!("{}", self.inner)
    }
}
