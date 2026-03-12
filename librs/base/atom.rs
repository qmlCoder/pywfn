use pyo3::prelude::*;

use rswfn;

use crate::base::Mole;

#[pyclass]
#[derive(Clone)]
pub struct Stm {
    pub inner: rswfn::base::Stm,
}

#[pymethods]
impl Stm {
    #[new]
    pub fn new(ex: [f64; 3], ey: [f64; 3], ez: [f64; 3]) -> PyResult<Self> {
        Ok(Self {
            inner: rswfn::base::Stm::new(ex, ey, ez),
        })
    }
}

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

    pub fn add(&mut self, atm: usize, xyz: [f64; 3], crg: f64, stm: Stm) {
        self.inner.add(atm, xyz, crg, stm.inner);
    }

    pub fn get(&self, idx: usize) -> PyResult<Atom> {
        let atom = self.inner.get(idx).clone();
        Ok(Atom { inner: atom })
    }

    pub fn len(&self) -> usize {
        self.inner.len()
    }

    pub fn dist(&self, i: usize, j: usize) -> f64 {
        self.inner.dist(i, j)
    }

    pub fn ato_nums(&self, form: &str, mole: &Mole) -> Vec<usize> {
        self.inner.ato_nums(form, &mole.inner)
    }
}
