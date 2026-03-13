use pyo3::prelude::*;

use rswfn;

use crate::base::Mole;

#[pyclass]
#[derive(Clone)]
pub struct Stm {
    pub core: rswfn::base::Stm,
}

#[pymethods]
impl Stm {
    #[new]
    pub fn new(ex: [f64; 3], ey: [f64; 3], ez: [f64; 3]) -> PyResult<Self> {
        Ok(Self {
            core: rswfn::base::Stm::new(ex, ey, ez),
        })
    }
}

#[pyclass]
pub struct Atom {
    core: rswfn::base::Atom,
}

#[pymethods]
impl Atom {
    pub fn rad(&self) -> f64 {
        self.core.rad
    }

    pub fn sym(&self) -> String {
        self.core.sym.clone()
    }
}

#[pyclass]
pub struct Atoms {
    core: rswfn::base::Atoms,
}

#[pymethods]
impl Atoms {
    #[new]
    pub fn new() -> PyResult<Self> {
        Ok(Self {
            core: rswfn::base::Atoms::new(),
        })
    }

    pub fn add(&mut self, atm: usize, xyz: [f64; 3], crg: f64, stm: Stm) {
        self.core.add(atm, xyz, crg, stm.core);
    }

    pub fn get(&self, idx: usize) -> PyResult<Atom> {
        let atom = self.core.get(idx).clone();
        Ok(Atom { core: atom })
    }

    pub fn len(&self) -> usize {
        self.core.len()
    }

    pub fn dist(&self, i: usize, j: usize) -> f64 {
        self.core.dist(i, j)
    }

    pub fn ato_nums(&self, form: &str, mole: &Mole) -> Vec<usize> {
        self.core.ato_nums(form, &mole.core)
    }
}
