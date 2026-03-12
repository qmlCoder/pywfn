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

    fn build(&mut self, data: Vec<(usize, usize, usize, f64, f64, (usize, usize, isize))>) {
        let mut basis_datas = vec![];
        for (atm, shl, ang, alp, coe, nlm) in data {
            basis_datas.push(rswfn::base::BasisData::new(atm, shl, ang, alp, coe, nlm));
        }
        self.inner.build(basis_datas);
    }

    pub fn __repr__(&self) -> String {
        format!("{}", self.inner)
    }
}
