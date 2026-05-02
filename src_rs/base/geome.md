```rust
use pyo3::prelude::*;

use rswfn;

use crate::base::Stm;

#[derive(Clone)]
#[pyclass]
pub struct Geome {
    pub core: rswfn::base::Geome,
}

#[pymethods]
impl Geome {
    #[new]
    pub fn new() -> PyResult<Self> {
        Ok(Self {
            core: rswfn::base::Geome::blank(),
        })
    }

    pub fn build(&mut self, atms: Vec<usize>, xyzs: Vec<[f64; 3]>, crgs: Vec<f64>, stms: Vec<Stm>) {
        let stms = stms.into_iter().map(|stm| stm.core).collect();
        self.core.build(&atms, &xyzs, &crgs, &stms);
    }

    pub fn __repr__(&self) -> String {
        format!("{}", self.core)
    }
}
```
