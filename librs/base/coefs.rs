use pyo3::prelude::*;

use numpy::PyReadonlyArray2;
use rswfn;

#[derive(Clone)]
#[pyclass]
pub struct Coefs {
    pub inner: rswfn::base::Coefs,
}

#[pymethods]
impl Coefs {
    #[new]
    fn new() -> Coefs {
        Coefs {
            inner: rswfn::base::Coefs::blank(),
        }
    }

    pub fn build(
        &mut self,
        ato_atms: Vec<usize>,
        ato_shls: Vec<usize>,
        ato_syms: Vec<String>,
        obt_engs: Vec<f64>,
        obt_occs: Vec<usize>,
        cmat: PyReadonlyArray2<f64>,
    ) {
        let cmat = cmat.as_array().to_owned();
        self.inner
            .build(ato_atms, ato_shls, ato_syms, obt_engs, obt_occs, cmat);
    }
    pub fn __repr__(&self) -> String {
        format!("{}", self.inner)
    }
}
