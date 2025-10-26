use pyo3::prelude::*;

use rswfn;

#[derive(Clone)]
#[pyclass]
pub struct Basis {
    pub inner: rswfn::base::Basis,
}

#[pymethods]
impl Basis {
    #[new]
    fn new() -> Self {
        Self {
            inner: rswfn::base::Basis::blank(),
        }
    }

    fn build(&mut self, data: Vec<(u32, u32, u32, f64, f64)>) {
        self.inner.build(data);
    }

    pub fn __repr__(&self) -> String {
        format!("{}", self.inner)
    }
}
