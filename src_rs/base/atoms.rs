use pyo3::prelude::*;

use rswfn;

use crate::base::Mole;

#[pyclass(from_py_object)]
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

    pub fn __repr__(&self) -> String {
        format!("{}", self.core)
    }
}

#[pyclass]
pub struct Atom {
    core: rswfn::base::Atom,
}

#[pymethods]
impl Atom {
    /// 添加一个原子
    pub fn rad(&self) -> f64 {
        self.core.rad
    }

    pub fn sym(&self) -> String {
        self.core.sym.clone()
    }
}

#[pyclass]
pub struct Atoms {
    pub core: rswfn::base::Atoms,
}

#[pymethods]
impl Atoms {
    #[new]
    pub fn new(
        atms: Vec<usize>,
        xyzs: Vec<[f64; 3]>,
        crgs: Vec<f64>,
        stms: Vec<Stm>,
    ) -> PyResult<Self> {
        let stms = match stms.len() {
            0 => None,
            _ => Some(stms.into_iter().map(|stm| stm.core).collect()),
        };
        Ok(Self {
            core: rswfn::base::Atoms::new(&atms, &xyzs, &crgs, stms.as_ref()),
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

    pub fn syms(&self) -> Vec<String> {
        self.core.syms()
    }

    pub fn xyzs(&self) -> Vec<[f64; 3]> {
        self.core.xyzs().to_vec()
    }

    pub fn __repr__(&self) -> String {
        format!("{}", self.core)
    }
}
