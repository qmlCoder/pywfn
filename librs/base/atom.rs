use pyo3::prelude::*;

use rswfn;

use crate::base::Mole;

#[pyclass]
pub struct Atom {
    inner: rswfn::base::Atom,
}

#[pymethods]
impl Atom {
    pub fn rad(&self) -> f64 {
        self.inner.rad
    }

    pub fn sym(&self) -> String {
        self.inner.sym.clone()
    }
}

#[pyclass]
pub struct Atoms {
    inner: rswfn::base::Atoms,
}

#[pymethods]
impl Atoms {
    #[new]
    pub fn new() -> PyResult<Self> {
        Ok(Self {
            inner: rswfn::base::Atoms::new(),
        })
    }

    pub fn add(&mut self, atm: u32, xyz: [f64; 3]) {
        self.inner.add(atm, xyz);
    }

    pub fn get(&self, idx: u32) -> PyResult<Atom> {
        let atom = self.inner.get(idx).clone();
        Ok(Atom { inner: atom })
    }

    pub fn len(&self) -> usize {
        self.inner.len()
    }

    pub fn dist(&self, i: u32, j: u32) -> f64 {
        self.inner.dist(i, j)
    }

    pub fn ato_nums(&self, form: &str, mole: &Mole) -> Vec<u32> {
        self.inner.ato_nums(form, &mole.inner)
    }
}
