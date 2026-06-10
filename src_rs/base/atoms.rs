use numpy::{IntoPyArray, PyArray2};
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

    pub fn set(&mut self, ex: [f64; 3], ey: [f64; 3], ez: [f64; 3]) {
        self.core.ex = ex;
        self.core.ey = ey;
        self.core.ez = ez;
    }

    #[getter]
    pub fn ex(&self) -> [f64; 3] {
        self.core.ex
    }

    #[getter]
    pub fn ey(&self) -> [f64; 3] {
        self.core.ey
    }

    #[getter]
    pub fn ez(&self) -> [f64; 3] {
        self.core.ez
    }

    #[getter]
    pub fn mat(&self, py: Python) -> Py<PyArray2<f64>> {
        self.core.mat.clone().into_pyarray(py).unbind()
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
    #[getter]
    pub fn rad(&self) -> f64 {
        self.core.rad
    }

    #[getter]
    pub fn sym(&self) -> String {
        self.core.sym.clone()
    }
}

#[derive(Clone)]
#[pyclass]
pub struct Atoms {
    pub core: rswfn::base::Atoms,
}

#[pymethods]
impl Atoms {
    #[new]
    pub fn new(atms: Vec<usize>, xyzs: Vec<[f64; 3]>, crgs: Vec<f64>) -> PyResult<Self> {
        Ok(Self {
            core: rswfn::base::Atoms::new(&atms, &xyzs, &crgs),
        })
    }

    pub fn add(&mut self, atm: usize, xyz: [f64; 3], crg: f64) {
        self.core.add(atm, xyz, crg);
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

    pub fn syms(&self) -> Vec<String> {
        self.core.syms()
    }

    pub fn xyzs(&self) -> Vec<[f64; 3]> {
        self.core.xyzs().to_vec()
    }

    pub fn nucs(&self) -> Vec<usize> {
        self.core.nucs().clone()
    }

    pub fn crgs(&self) -> Vec<f64> {
        self.core.crgs().clone()
    }

    pub fn nebs(&self, atm: usize) -> Vec<usize> {
        self.core.nebs(atm).clone()
    }

    pub fn is_line(&self, atm: usize) -> bool {
        self.core.is_line(atm)
    }

    pub fn is_bond(&self, i: usize, j: usize) -> bool {
        self.core.is_bond(i, j)
    }

    pub fn __repr__(&self) -> String {
        format!("{}", self.core)
    }
}
